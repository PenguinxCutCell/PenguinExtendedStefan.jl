function l2_error(values_num::AbstractVector{T}, values_ref::AbstractVector{T}, weights::AbstractVector{T}) where {T<:Real}
    n = length(values_num)
    length(values_ref) == n || throw(DimensionMismatch("values_ref length mismatch"))
    length(weights) == n || throw(DimensionMismatch("weights length mismatch"))

    acc = zero(T)
    wsum = zero(T)
    @inbounds for i in 1:n
        w = weights[i]
        if !(isfinite(w) && w > zero(T))
            continue
        end
        d = values_num[i] - values_ref[i]
        acc += w * d^2
        wsum += w
    end
    return (error=sqrt(acc), weight=wsum)
end

function linf_error(values_num::AbstractVector{T}, values_ref::AbstractVector{T}) where {T<:Real}
    n = length(values_num)
    length(values_ref) == n || throw(DimensionMismatch("values_ref length mismatch"))
    m = zero(T)
    @inbounds for i in 1:n
        m = max(m, abs(values_num[i] - values_ref[i]))
    end
    return m
end

function profile_error_on_phase(
    values_num::AbstractVector{T},
    values_ref,
    cap::AssembledCapacity{N,T},
    t::T;
    tol::T=zero(T),
    use_interface_band::Bool=false,
) where {N,T}
    length(values_num) == cap.ntotal || throw(DimensionMismatch("values_num length must be $(cap.ntotal)"))

    ref = Vector{T}(undef, cap.ntotal)
    weights = zeros(T, cap.ntotal)

    @inbounds for i in 1:cap.ntotal
        w = cap.buf.V[i]
        if !(isfinite(w) && w > tol)
            ref[i] = zero(T)
            continue
        end
        if use_interface_band && !(isfinite(cap.buf.Γ[i]) && cap.buf.Γ[i] > tol)
            ref[i] = zero(T)
            continue
        end
        ref[i] = _eval_exact(values_ref, cap.C_ω[i], t)
        weights[i] = w
    end

    l2 = l2_error(values_num, ref, weights)
    linf = linf_error(values_num, ref)
    return (L2=l2.error, Linf=linf, weight=l2.weight)
end

function interface_position_from_levelset_1d(rep::AbstractInterfaceRep{1,T}, grid::CartesianGrids.CartesianGrid{1,T}) where {T}
    x = CartesianGrids.grid1d(grid, 1)
    ϕ = vec(phi_values(rep))
    for i in 1:(length(x) - 1)
        a = ϕ[i]
        b = ϕ[i + 1]
        if a == zero(T)
            return x[i]
        elseif a * b < zero(T)
            θ = abs(a) / (abs(a) + abs(b))
            return (one(T) - θ) * x[i] + θ * x[i + 1]
        end
    end
    return T(NaN)
end

interface_position_error(num, ref) = abs(num - ref)

stefan_residual_metrics(solver::ExtendedStefanSolver) = extended_stefan_residual_metrics(solver)
