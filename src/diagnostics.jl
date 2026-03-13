function active_mask(cap::AssembledCapacity{N,T}; tol::T=zero(T)) where {N,T}
    m = falses(cap.ntotal)
    @inbounds for i in eachindex(m)
        v = cap.buf.V[i]
        m[i] = isfinite(v) && v > tol
    end
    return m
end

function phase_volume(cap::AssembledCapacity{N,T}; tol::T=zero(T)) where {N,T}
    s = zero(T)
    @inbounds for i in 1:cap.ntotal
        v = cap.buf.V[i]
        if isfinite(v) && v > tol
            s += v
        end
    end
    return s
end

function phase_volume(V::AbstractVector{T}; tol::T=zero(T)) where {T<:Real}
    s = zero(T)
    @inbounds for i in eachindex(V)
        v = V[i]
        if isfinite(v) && v > tol
            s += v
        end
    end
    return s
end

function _eval_exact(Texact, x::SVector{N,T}, t::T) where {N,T}
    if applicable(Texact, x..., t)
        return convert(T, Texact(x..., t))
    elseif applicable(Texact, x...)
        return convert(T, Texact(x...))
    end
    throw(ArgumentError("exact solution callback must accept (x...) or (x..., t)"))
end

function _component_error_norms(
    uω::AbstractVector{T},
    uexact,
    cap::AssembledCapacity{N,T},
    t::T;
    tol::T=zero(T),
) where {N,T}
    length(uω) == cap.ntotal || throw(DimensionMismatch("uω length must be $(cap.ntotal)"))

    l1 = zero(T)
    l2 = zero(T)
    linf = zero(T)
    wsum = zero(T)

    @inbounds for i in 1:cap.ntotal
        w = cap.buf.V[i]
        if !(isfinite(w) && w > tol)
            continue
        end
        e = abs(uω[i] - _eval_exact(uexact, cap.C_ω[i], t))
        l1 += w * e
        l2 += w * e^2
        linf = max(linf, e)
        wsum += w
    end

    return (L1=l1, L2=sqrt(l2), Linf=linf, weight=wsum)
end

function temperature_error_norms(
    Tlω::AbstractVector{T},
    Tsω::AbstractVector{T},
    Texact_l,
    Texact_s,
    cap_l::AssembledCapacity{N,T},
    cap_s::AssembledCapacity{N,T},
    t::T;
    tol::T=zero(T),
) where {N,T}
    e_l = _component_error_norms(Tlω, Texact_l, cap_l, t; tol=tol)
    e_s = _component_error_norms(Tsω, Texact_s, cap_s, t; tol=tol)
    return (liquid=e_l, solid=e_s)
end

function concentration_error_norms(
    Clω::AbstractVector{T},
    Csω::AbstractVector{T},
    Cexact_l,
    Cexact_s,
    cap_l::AssembledCapacity{N,T},
    cap_s::AssembledCapacity{N,T},
    t::T;
    tol::T=zero(T),
) where {N,T}
    e_l = _component_error_norms(Clω, Cexact_l, cap_l, t; tol=tol)
    e_s = _component_error_norms(Csω, Cexact_s, cap_s, t; tol=tol)
    return (liquid=e_l, solid=e_s)
end

function extended_stefan_residual_metrics(solver::ExtendedStefanSolver)
    cache = solver.cache
    state = solver.state
    return (
        thermal_stefan_residual_norm=cache.metrics[:stefan_residual_norm],
        solute_residual_norm=cache.metrics[:solute_residual_norm],
        speed_increment_norm=cache.metrics[:speed_increment_norm],
        max_speed=cache.metrics[:max_speed],
        nonlinear_iteration_count=cache.metrics[:nonlinear_iters],
        last_time=state.t,
    )
end
