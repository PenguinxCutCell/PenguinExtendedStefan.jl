function _eval_state_value(v, x::SVector{N,T}, t::T) where {N,T}
    if v isa Number
        return convert(T, v)
    elseif v isa Function
        if applicable(v, x..., t)
            return convert(T, v(x..., t))
        elseif applicable(v, x...)
            return convert(T, v(x...))
        end
    end
    throw(ArgumentError("state initializer must be scalar, vector, or callback (x...) / (x..., t)"))
end

function _init_component_values(v0, xpts::Vector{SVector{N,T}}, t::T, default_v) where {N,T}
    nt = length(xpts)
    if v0 === nothing
        out = Vector{T}(undef, nt)
        @inbounds for i in 1:nt
            out[i] = _eval_state_value(default_v, xpts[i], t)
        end
        return out
    elseif v0 isa AbstractVector
        length(v0) == nt || throw(DimensionMismatch("initializer length must be $nt"))
        return convert(Vector{T}, v0)
    end

    out = Vector{T}(undef, nt)
    @inbounds for i in 1:nt
        out[i] = _eval_state_value(v0, xpts[i], t)
    end
    return out
end

function _init_speed(speed0, ::Type{T}, shp::NTuple{N,Int}) where {N,T}
    if speed0 === nothing
        return zeros(T, shp)
    elseif speed0 isa Number
        return fill(convert(T, speed0), shp)
    elseif speed0 isa AbstractArray
        size(speed0) == shp || throw(DimensionMismatch("speed0 shape must be $(shp)"))
        return convert(Array{T,N}, speed0)
    end
    throw(ArgumentError("speed0 must be nothing, scalar, or array with shape grid.n"))
end

function build_solver(
    prob::ExtendedStefanProblem{N,T};
    t0::T=zero(T),
    Tlω0=nothing,
    Tlγ0=nothing,
    Tsω0=nothing,
    Tsγ0=nothing,
    Clω0=nothing,
    Clγ0=nothing,
    Csω0=nothing,
    Csγ0=nothing,
    speed0=nothing,
) where {N,T}
    cache = build_cache(prob)

    ϕ0 = phi_values(prob.interface_rep)
    update_slab_field!(cache.slab, ϕ0, ϕ0, t0, t0 + one(T))
    PenguinDiffusion._build_moving_slab!(cache.modelT, t0, one(T))
    PenguinDiffusion._build_moving_slab!(cache.modelC, t0, one(T))

    cap_l = something(cache.modelT.cap1_slab)
    cap_s = something(cache.modelT.cap2_slab)

    Tm_default = prob.params.alloy_ic.T_m

    Tlω = _init_component_values(Tlω0, cap_l.C_ω, t0, Tm_default)
    Tlγ = _init_component_values(Tlγ0, cap_l.C_γ, t0, Tm_default)
    Tsω = _init_component_values(Tsω0, cap_s.C_ω, t0, Tm_default)
    Tsγ = _init_component_values(Tsγ0, cap_s.C_γ, t0, Tm_default)

    zero_cb = (args...) -> zero(T)
    Clω = _init_component_values(Clω0, cap_l.C_ω, t0, zero_cb)
    Clγ = _init_component_values(Clγ0, cap_l.C_γ, t0, zero_cb)
    Csω = _init_component_values(Csω0, cap_s.C_ω, t0, zero_cb)
    Csγ = _init_component_values(Csγ0, cap_s.C_γ, t0, zero_cb)

    speed = _init_speed(speed0, T, size(ϕ0))
    frozen = falses(size(ϕ0))

    logs = Dict{Symbol,Any}(
        :step => 0,
        :times => T[t0],
        :speed_max => T[],
        :nonlinear_iters => Int[],
        :solute_residual => T[],
        :thermal_residual => T[],
        :speed_increment => T[],
    )

    state = ExtendedStefanState{N,T,typeof(speed)}(
        t0,
        Tlω,
        Tlγ,
        Tsω,
        Tsγ,
        Clω,
        Clγ,
        Csω,
        Csγ,
        speed,
        frozen,
        logs,
    )

    return ExtendedStefanSolver{N,T,typeof(prob),typeof(state),typeof(cache)}(prob, state, cache)
end

function _convergence_reached(
    speed_inc,
    sol_res,
    th_res,
    vref,
    opts::ExtendedStefanOptions{T},
) where {T}
    atol = opts.coupling_tol
    rtol = opts.coupling_reltol
    thresh = atol + rtol * max(one(T), vref)
    return speed_inc <= thresh && sol_res <= thresh && th_res <= thresh
end

function step!(
    solver::ExtendedStefanSolver{N,T},
    dt::T;
    method::Symbol=:linearsolve,
    kwargs...,
) where {N,T}
    dt > zero(T) || throw(ArgumentError("dt must be positive"))

    prob = solver.problem
    state = solver.state
    cache = solver.cache
    rep = prob.interface_rep

    t = state.t
    phi_n = phi_values(rep)

    Tlω_prev = copy(state.Tlω)
    Tlγ_prev = copy(state.Tlγ)
    Tsω_prev = copy(state.Tsω)
    Tsγ_prev = copy(state.Tsγ)
    Clω_prev = copy(state.Clω)
    Clγ_prev = copy(state.Clγ)
    Csω_prev = copy(state.Csω)
    Csγ_prev = copy(state.Csγ)

    Vlag = copy(state.speed_full)

    Tlω_new = copy(state.Tlω)
    Tlγ_new = copy(state.Tlγ)
    Tsω_new = copy(state.Tsω)
    Tsγ_new = copy(state.Tsγ)
    Clω_new = copy(state.Clω)
    Clγ_new = copy(state.Clγ)
    Csω_new = copy(state.Csω)
    Csγ_new = copy(state.Csγ)

    Vnew = similar(Vlag)
    frozen_new = similar(state.frozen_mask)

    iter_done = 0
    speed_inc = typemax(T)
    sol_res_norm = typemax(T)
    th_res_norm = typemax(T)

    for iter in 1:prob.options.coupling_max_iter
        iter_done = iter

        phi_pred = predict_phi(rep, Vlag, t, dt)

        Tlω_new,
        Tlγ_new,
        Tsω_new,
        Tsγ_new,
        Clω_new,
        Clγ_new,
        Csω_new,
        Csγ_new,
        _, _ = solve_monolithic_system!(
            cache,
            phi_n,
            phi_pred,
            Vlag,
            Tlω_prev,
            Tlγ_prev,
            Tsω_prev,
            Tsγ_prev,
            Clω_prev,
            Clγ_prev,
            Csω_prev,
            Csγ_prev,
            t,
            dt,
            prob;
            method=method,
            kwargs...,
        )

        stefan_speed!(
            Vnew,
            frozen_new,
            cache,
            Tlω_new,
            Tlγ_new,
            Tsω_new,
            Tsγ_new,
            Tlω_prev,
            Tlγ_prev,
            Tsω_prev,
            Tsγ_prev,
            prob.params.rhoL,
            prob.options.scheme,
            t,
            dt,
        )

        sol_res = solute_conservation_residual(
            cache,
            Clω_new,
            Clγ_new,
            Csω_new,
            Csγ_new,
            Vnew,
            prob.options.scheme,
            t,
            dt,
        )
        th_res = thermal_stefan_residual(
            cache,
            Tlω_new,
            Tlγ_new,
            Tsω_new,
            Tsγ_new,
            Tlω_prev,
            Tlγ_prev,
            Tsω_prev,
            Tsγ_prev,
            Vnew,
            prob.params.rhoL,
            prob.options.scheme,
            t,
            dt,
        )

        speed_inc = maximum(abs.(vec(Vnew .- Vlag)))
        sol_res_norm = sol_res.maxabs
        th_res_norm = th_res.maxabs

        ω = convert(T, prob.options.coupling_damping)
        Vlag .= (one(T) - ω) .* Vlag .+ ω .* Vnew

        vref = maximum(abs.(vec(Vlag)))
        if _convergence_reached(speed_inc, sol_res_norm, th_res_norm, vref, prob.options)
            break
        end
    end

    extend_speed!(
        rep,
        Vlag;
        frozen=frozen_new,
        nb_iters=prob.options.extend_iters,
        cfl=prob.options.extend_cfl,
        interface_band=prob.options.interface_band,
    )
    advance!(rep, Vlag, t, dt)

    step_id = get(state.logs, :step, 0) + 1
    state.logs[:step] = step_id
    if prob.options.reinit && prob.options.reinit_every > 0 && (step_id % prob.options.reinit_every == 0)
        reinit!(rep)
    end

    state.t = t + dt
    state.Tlω .= Tlω_new
    state.Tlγ .= Tlγ_new
    state.Tsω .= Tsω_new
    state.Tsγ .= Tsγ_new
    state.Clω .= Clω_new
    state.Clγ .= Clγ_new
    state.Csω .= Csω_new
    state.Csγ .= Csγ_new
    state.speed_full .= Vlag
    state.frozen_mask .= frozen_new

    max_speed = maximum(abs.(vec(Vlag)))
    push!(state.logs[:times], state.t)
    push!(state.logs[:speed_max], max_speed)
    push!(state.logs[:nonlinear_iters], iter_done)
    push!(state.logs[:solute_residual], sol_res_norm)
    push!(state.logs[:thermal_residual], th_res_norm)
    push!(state.logs[:speed_increment], speed_inc)

    cache.metrics[:solute_residual_norm] = sol_res_norm
    cache.metrics[:stefan_residual_norm] = th_res_norm
    cache.metrics[:speed_increment_norm] = speed_inc
    cache.metrics[:max_speed] = max_speed
    cache.metrics[:nonlinear_iters] = iter_done

    return state
end

function solve!(
    solver::ExtendedStefanSolver{N,T},
    tspan::Tuple{T,T};
    dt::T,
    save_history::Bool=true,
    method::Symbol=:linearsolve,
    kwargs...,
) where {N,T}
    t0, tend = tspan
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))
    dt > zero(T) || throw(ArgumentError("dt must be positive"))

    state = solver.state
    if abs(state.t - t0) > sqrt(eps(T)) * max(one(T), abs(t0))
        throw(ArgumentError("solver state time $(state.t) must match tspan start $t0"))
    end

    history = Vector{NamedTuple}()
    if save_history
        push!(history, (
            t=state.t,
            Tlω=copy(state.Tlω), Tlγ=copy(state.Tlγ),
            Tsω=copy(state.Tsω), Tsγ=copy(state.Tsγ),
            Clω=copy(state.Clω), Clγ=copy(state.Clγ),
            Csω=copy(state.Csω), Csγ=copy(state.Csγ),
            speed=copy(state.speed_full),
            phi=phi_values(solver.problem.interface_rep),
        ))
    end

    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))
    while state.t < tend - tol
        dt_step = min(dt, tend - state.t)
        step!(solver, dt_step; method=method, kwargs...)
        if save_history
            push!(history, (
                t=state.t,
                Tlω=copy(state.Tlω), Tlγ=copy(state.Tlγ),
                Tsω=copy(state.Tsω), Tsγ=copy(state.Tsγ),
                Clω=copy(state.Clω), Clγ=copy(state.Clγ),
                Csω=copy(state.Csω), Csγ=copy(state.Csγ),
                speed=copy(state.speed_full),
                phi=phi_values(solver.problem.interface_rep),
            ))
        end
    end

    return (
        solver=solver,
        state=solver.state,
        times=copy(solver.state.logs[:times]),
        history=history,
    )
end
