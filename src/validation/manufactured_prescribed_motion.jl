function prescribed_planar_motion_case(
    ;
    h0::Real=0.35,
    V0::Real=0.15,
    ω::Real=0.0,
    A::Real=0.0,
    k_partition::Real=0.7,
    Tm::Real=0.0,
    m_liquidus::Real=-0.2,
    rho_cp_l::Real=1.0,
    rho_cp_s::Real=1.0,
    kappa_l::Real=1.0,
    kappa_s::Real=0.8,
    D_l::Real=0.1,
    D_s::Real=0.04,
)
    k = float(k_partition)
    Tm0 = float(Tm)
    mliq = float(m_liquidus)
    h(t) = float(h0) + float(V0) * t + float(A) * sin(float(ω) * t)
    V(t) = float(V0) + float(A) * float(ω) * cos(float(ω) * t)

    ClΓ(t) = 0.22 + 0.03 * sin(t)
    dClΓdt(t) = 0.03 * cos(t)
    CsΓ(t) = k * ClΓ(t)
    TΓ(t) = Tm0 + mliq * ClΓ(t)
    dTΓdt(t) = mliq * dClΓdt(t)

    ξ(x, t) = x - h(t)

    Tl_exact(x, t) = TΓ(t) + 0.20 * ξ(x, t) + 0.07 * ξ(x, t)^2
    Ts_exact(x, t) = TΓ(t) - 0.14 * ξ(x, t) + 0.05 * ξ(x, t)^2
    Cl_exact(x, t) = ClΓ(t) + 0.09 * ξ(x, t) + 0.04 * ξ(x, t)^2
    Cs_exact(x, t) = CsΓ(t) - 0.07 * ξ(x, t) + 0.03 * ξ(x, t)^2

    dξdx = 1.0
    source_T_l(x, t) = rho_cp_l * (dTΓdt(t) - V(t) * (0.20 * dξdx + 2 * 0.07 * ξ(x, t))) - kappa_l * (2 * 0.07)
    source_T_s(x, t) = rho_cp_s * (dTΓdt(t) - V(t) * (-0.14 * dξdx + 2 * 0.05 * ξ(x, t))) - kappa_s * (2 * 0.05)
    source_C_l(x, t) = (dClΓdt(t) - V(t) * (0.09 * dξdx + 2 * 0.04 * ξ(x, t))) - D_l * (2 * 0.04)
    source_C_s(x, t) = (k * dClΓdt(t) - V(t) * (-0.07 * dξdx + 2 * 0.03 * ξ(x, t))) - D_s * (2 * 0.03)

    return (
        h=h,
        V=V,
        phi0=(x) -> x - h(0.0),
        speed=(x, t) -> V(t),
        Tl_exact=Tl_exact,
        Ts_exact=Ts_exact,
        Cl_exact=Cl_exact,
        Cs_exact=Cs_exact,
        source_T_l=source_T_l,
        source_T_s=source_T_s,
        source_C_l=source_C_l,
        source_C_s=source_C_s,
    )
end

function prescribed_circle_motion_case(
    ;
    center::Tuple{Real,Real}=(0.0, 0.0),
    R0::Real=0.30,
    V0::Real=0.05,
    k_partition::Real=0.75,
    Tm::Real=0.0,
    m_liquidus::Real=-0.2,
    rho_cp_l::Real=1.0,
    rho_cp_s::Real=1.0,
    kappa_l::Real=1.0,
    kappa_s::Real=0.9,
    D_l::Real=0.1,
    D_s::Real=0.04,
)
    xc, yc = float(center[1]), float(center[2])
    k = float(k_partition)
    Tm0 = float(Tm)
    mliq = float(m_liquidus)

    R(t) = float(R0) + float(V0) * t
    V(t) = float(V0)

    δ(x, y, t) = hypot(x - xc, y - yc) - R(t)
    ξ(x) = (x - xc) / max(R0, eps())

    ClΓ(x, y, t) = 0.20 + 0.02 * sin(t) + 0.01 * ξ(x)
    CsΓ(x, y, t) = k * ClΓ(x, y, t)
    TΓ(x, y, t) = Tm0 + mliq * ClΓ(x, y, t)

    Tl_exact(x, y, t) = TΓ(x, y, t) + 0.15 * δ(x, y, t) + 0.05 * δ(x, y, t)^2
    Ts_exact(x, y, t) = TΓ(x, y, t) - 0.10 * δ(x, y, t) + 0.04 * δ(x, y, t)^2
    Cl_exact(x, y, t) = ClΓ(x, y, t) + 0.08 * δ(x, y, t) + 0.03 * δ(x, y, t)^2
    Cs_exact(x, y, t) = CsΓ(x, y, t) - 0.06 * δ(x, y, t) + 0.02 * δ(x, y, t)^2

    source_T_l(x, y, t) = rho_cp_l * _fd_dt_2d(Tl_exact, x, y, t) - kappa_l * _fd_lap_2d(Tl_exact, x, y, t)
    source_T_s(x, y, t) = rho_cp_s * _fd_dt_2d(Ts_exact, x, y, t) - kappa_s * _fd_lap_2d(Ts_exact, x, y, t)
    source_C_l(x, y, t) = _fd_dt_2d(Cl_exact, x, y, t) - D_l * _fd_lap_2d(Cl_exact, x, y, t)
    source_C_s(x, y, t) = _fd_dt_2d(Cs_exact, x, y, t) - D_s * _fd_lap_2d(Cs_exact, x, y, t)

    return (
        R=R,
        V=V,
        phi0=(x, y) -> hypot(x - xc, y - yc) - R(0.0),
        speed=(x, y, t) -> V(t),
        Tl_exact=Tl_exact,
        Ts_exact=Ts_exact,
        Cl_exact=Cl_exact,
        Cs_exact=Cs_exact,
        source_T_l=source_T_l,
        source_T_s=source_T_s,
        source_C_l=source_C_l,
        source_C_s=source_C_s,
    )
end

function _speed_nodes(grid::CartesianGrids.CartesianGrid{N,T}, speed_fun, t::T) where {N,T}
    xyz = CartesianGrids.grid1d(grid)
    out = Array{T,N}(undef, grid.n...)
    @inbounds for I in CartesianIndices(out)
        x = ntuple(d -> convert(T, xyz[d][I[d]]), N)
        out[I] = if applicable(speed_fun, x..., t)
            convert(T, speed_fun(x..., t))
        elseif applicable(speed_fun, x...)
            convert(T, speed_fun(x...))
        else
            throw(ArgumentError("speed callback must accept (x..., t) or (x...)"))
        end
    end
    return out
end

function _frozen_from_capacity(cap::AssembledCapacity{N,T}; tol::T=sqrt(eps(T))) where {N,T}
    mask = falses(cap.nnodes)
    mflat = vec(mask)
    Γ = cap.buf.Γ
    @inbounds for i in eachindex(Γ)
        γ = Γ[i]
        mflat[i] = isfinite(γ) && γ > tol
    end
    return mask
end

function solve_prescribed_motion!(
    solver::ExtendedStefanSolver{N,T},
    tspan::Tuple{T,T};
    dt::T,
    speed_callback,
    save_history::Bool=true,
    method::Symbol=:linearsolve,
    kwargs...,
) where {N,T}
    t0, tend = tspan
    state = solver.state
    prob = solver.problem
    cache = solver.cache
    rep = prob.interface_rep

    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    if abs(state.t - t0) > sqrt(eps(T)) * max(one(T), abs(t0))
        throw(ArgumentError("solver state time $(state.t) must match tspan start $t0"))
    end

    history = Vector{NamedTuple}()
    if save_history
        push!(history, (t=state.t, phi=phi_values(rep), speed=copy(state.speed_full)))
    end

    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))
    while state.t < tend - tol
        dt_step = min(dt, tend - state.t)
        t = state.t

        phi_n = phi_values(rep)
        Vprescribed = _speed_nodes(prob.grid, speed_callback, t)
        phi_pred = predict_phi(rep, Vprescribed, t, dt_step)

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
            Vprescribed,
            state.Tlω,
            state.Tlγ,
            state.Tsω,
            state.Tsγ,
            state.Clω,
            state.Clγ,
            state.Csω,
            state.Csγ,
            t,
            dt_step,
            prob;
            method=method,
            kwargs...,
        )

        cap = something(cache.modelT.cap1_slab)
        frozen = _frozen_from_capacity(cap)

        extend_speed!(
            rep,
            Vprescribed;
            frozen=frozen,
            nb_iters=prob.options.extend_iters,
            cfl=prob.options.extend_cfl,
            interface_band=prob.options.interface_band,
        )
        advance!(rep, Vprescribed, t, dt_step)

        state.t = t + dt_step
        state.Tlω .= Tlω_new
        state.Tlγ .= Tlγ_new
        state.Tsω .= Tsω_new
        state.Tsγ .= Tsγ_new
        state.Clω .= Clω_new
        state.Clγ .= Clγ_new
        state.Csω .= Csω_new
        state.Csγ .= Csγ_new
        state.speed_full .= Vprescribed
        state.frozen_mask .= frozen

        if save_history
            push!(history, (t=state.t, phi=phi_values(rep), speed=copy(state.speed_full)))
        end
    end

    return (solver=solver, state=solver.state, history=history)
end
