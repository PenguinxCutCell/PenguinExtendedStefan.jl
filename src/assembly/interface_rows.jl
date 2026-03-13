@inline function _finite_point(x::SVector)
    return all(isfinite, x)
end

function _point_values(v, xpts::AbstractVector{<:SVector{N,T}}, t::T) where {N,T}
    out = Vector{T}(undef, length(xpts))
    @inbounds for i in eachindex(xpts)
        x = xpts[i]
        out[i] = _finite_point(x) ? convert(T, eval_bc(v, x, t)) : zero(T)
    end
    return out
end

function _alloy_trace_values(prob::ExtendedStefanProblem{N,T}, Cγ, t::T) where {N,T}
    ic = prob.params.alloy_ic
    k = _point_values(ic.k_partition, Cγ, t)
    Tm = _point_values(ic.T_m, Cγ, t)
    m = _point_values(ic.m_liquidus, Cγ, t)
    return k, Tm, m
end

function _active_rows_extended(
    cap_l::AssembledCapacity{N,T},
    cap_s::AssembledCapacity{N,T},
    lay::ExtendedLayout;
    tol::T=zero(T),
) where {N,T}
    activeω_l, activeγ_l = PenguinDiffusion._cell_activity_masks(cap_l)
    activeω_s, activeγ_s = PenguinDiffusion._cell_activity_masks(cap_s)

    nt = lay.nt
    active = falses(8 * nt)
    @inbounds for i in 1:nt
        active[lay.Tlω[i]] = activeω_l[i]
        active[lay.Tsω[i]] = activeω_s[i]
        active[lay.Clω[i]] = activeω_l[i]
        active[lay.Csω[i]] = activeω_s[i]

        γok = activeγ_l[i] && activeγ_s[i]
        active[lay.Tlγ[i]] = γok
        active[lay.Tsγ[i]] = γok
        active[lay.Clγ[i]] = γok
        active[lay.Csγ[i]] = γok
    end
    return active
end

function apply_interface_rows!(
    A::SparseMatrixCSC{T,Int},
    b::Vector{T},
    cache::ExtendedStefanCache{N,T},
    prob::ExtendedStefanProblem{N,T},
    cap_l::AssembledCapacity{N,T},
    cap_s::AssembledCapacity{N,T},
    JCl::SparseMatrixCSC{T,Int},
    LCl::SparseMatrixCSC{T,Int},
    JCs::SparseMatrixCSC{T,Int},
    LCs::SparseMatrixCSC{T,Int},
    Vlag_full,
    ttrace::T,
) where {N,T}
    lay = cache.layout
    nt = lay.nt

    k_part, Tm, mliq = _alloy_trace_values(prob, cap_l.C_γ, ttrace)
    vlag = vec(Vlag_full)
    length(vlag) == nt || throw(DimensionMismatch("Vlag shape must match grid.n (nt=$nt)"))

    @inbounds for i in 1:nt
        row_cont = lay.Tlγ[i]
        A[row_cont, lay.Tlγ[i]] = one(T)
        A[row_cont, lay.Tsγ[i]] = -one(T)
        b[row_cont] = zero(T)

        row_liq = lay.Tsγ[i]
        A[row_liq, lay.Tlγ[i]] = one(T)
        A[row_liq, lay.Clγ[i]] = -mliq[i]
        b[row_liq] = Tm[i]

        row_part = lay.Clγ[i]
        A[row_part, lay.Csγ[i]] = one(T)
        A[row_part, lay.Clγ[i]] = -k_part[i]
        b[row_part] = zero(T)

        row_sol = lay.Csγ[i]
        b[row_sol] = zero(T)

        for j in 1:nt
            v1 = JCs[i, j]
            iszero(v1) || (A[row_sol, lay.Csω[j]] = A[row_sol, lay.Csω[j]] + v1)
            v2 = LCs[i, j]
            iszero(v2) || (A[row_sol, lay.Csγ[j]] = A[row_sol, lay.Csγ[j]] + v2)

            v3 = JCl[i, j]
            iszero(v3) || (A[row_sol, lay.Clω[j]] = A[row_sol, lay.Clω[j]] - v3)
            v4 = LCl[i, j]
            iszero(v4) || (A[row_sol, lay.Clγ[j]] = A[row_sol, lay.Clγ[j]] - v4)
        end

        vi = convert(T, vlag[i])
        A[row_sol, lay.Clγ[i]] = A[row_sol, lay.Clγ[i]] - vi
        A[row_sol, lay.Csγ[i]] = A[row_sol, lay.Csγ[i]] + vi
    end

    return A, b
end
