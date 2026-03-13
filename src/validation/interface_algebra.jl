function _residual_norms(r::AbstractVector{T}) where {T<:Real}
    if isempty(r)
        return (maxabs=zero(T), l2=zero(T), count=0)
    end
    maxr = maximum(abs.(r))
    l2 = sqrt(sum(abs2, r) / length(r))
    return (maxabs=maxr, l2=l2, count=length(r))
end

function interface_algebra_residuals(
    state::ExtendedStefanState{N,T},
    params::ExtendedStefanParams;
    mask::Union{Nothing,AbstractVector{Bool}}=nothing,
) where {N,T}
    nt = length(state.Tlγ)
    if !(mask === nothing)
        length(mask) == nt || throw(DimensionMismatch("mask length must be $nt"))
    end

    k = params.alloy_ic.k_partition
    Tm = params.alloy_ic.T_m
    m = params.alloy_ic.m_liquidus

    r_cont = T[]
    r_part = T[]
    r_liq = T[]

    @inbounds for i in 1:nt
        if !(mask === nothing) && !mask[i]
            continue
        end
        push!(r_cont, state.Tlγ[i] - state.Tsγ[i])
        push!(r_part, state.Csγ[i] - k * state.Clγ[i])
        push!(r_liq, state.Tlγ[i] - (Tm + m * state.Clγ[i]))
    end

    return (
        continuity=_residual_norms(r_cont),
        partition=_residual_norms(r_part),
        liquidus=_residual_norms(r_liq),
        residual_continuity=r_cont,
        residual_partition=r_part,
        residual_liquidus=r_liq,
    )
end

function interface_algebra_residuals(
    solver::ExtendedStefanSolver{N,T};
    tol::T=sqrt(eps(T)),
) where {N,T}
    cap = solver.cache.modelC.cap1_slab
    cap === nothing && throw(ArgumentError("slab geometry not initialized"))
    Γ = cap.buf.Γ
    mask = BitVector(undef, length(Γ))
    @inbounds for i in eachindex(Γ)
        γ = Γ[i]
        mask[i] = isfinite(γ) && γ > tol
    end
    return interface_algebra_residuals(solver.state, solver.problem.params; mask=mask)
end
