mutable struct SlabPhiField{N,T,AT}
    xyz::NTuple{N,Vector{T}}
    phi_n::AT
    phi_np1::AT
    t0::T
    t1::T
end

@inline function _lerp(a::T, b::T, θ::T) where {T}
    return (one(T) - θ) * a + θ * b
end

function _interp1(xnodes::AbstractVector{T}, vals::AbstractVector{T}, x::T) where {T}
    n = length(xnodes)
    n == 1 && return vals[1]
    if x <= xnodes[1]
        return vals[1]
    elseif x >= xnodes[end]
        return vals[end]
    end
    i = clamp(searchsortedlast(xnodes, x), 1, n - 1)
    x0 = xnodes[i]
    x1 = xnodes[i + 1]
    θ = (x - x0) / (x1 - x0)
    return _lerp(vals[i], vals[i + 1], θ)
end

function _interp2(
    xnodes::AbstractVector{T},
    ynodes::AbstractVector{T},
    vals::AbstractMatrix{T},
    x::T,
    y::T,
) where {T}
    nx = length(xnodes)
    ny = length(ynodes)
    nx == 1 && ny == 1 && return vals[1, 1]

    i = x <= xnodes[1] ? 1 : (x >= xnodes[end] ? nx - 1 : clamp(searchsortedlast(xnodes, x), 1, nx - 1))
    j = y <= ynodes[1] ? 1 : (y >= ynodes[end] ? ny - 1 : clamp(searchsortedlast(ynodes, y), 1, ny - 1))
    x0, x1 = xnodes[i], xnodes[i + 1]
    y0, y1 = ynodes[j], ynodes[j + 1]
    θx = x1 == x0 ? zero(T) : (x - x0) / (x1 - x0)
    θy = y1 == y0 ? zero(T) : (y - y0) / (y1 - y0)

    v00 = vals[i, j]
    v10 = vals[i + 1, j]
    v01 = vals[i, j + 1]
    v11 = vals[i + 1, j + 1]

    vx0 = _lerp(v00, v10, θx)
    vx1 = _lerp(v01, v11, θx)
    return _lerp(vx0, vx1, θy)
end

function _interp3(
    xnodes::AbstractVector{T},
    ynodes::AbstractVector{T},
    znodes::AbstractVector{T},
    vals::Array{T,3},
    x::T,
    y::T,
    z::T,
) where {T}
    nx = length(xnodes)
    ny = length(ynodes)
    nz = length(znodes)
    nx == 1 && ny == 1 && nz == 1 && return vals[1, 1, 1]

    i = x <= xnodes[1] ? 1 : (x >= xnodes[end] ? nx - 1 : clamp(searchsortedlast(xnodes, x), 1, nx - 1))
    j = y <= ynodes[1] ? 1 : (y >= ynodes[end] ? ny - 1 : clamp(searchsortedlast(ynodes, y), 1, ny - 1))
    k = z <= znodes[1] ? 1 : (z >= znodes[end] ? nz - 1 : clamp(searchsortedlast(znodes, z), 1, nz - 1))

    x0, x1 = xnodes[i], xnodes[i + 1]
    y0, y1 = ynodes[j], ynodes[j + 1]
    z0, z1 = znodes[k], znodes[k + 1]
    θx = x1 == x0 ? zero(T) : (x - x0) / (x1 - x0)
    θy = y1 == y0 ? zero(T) : (y - y0) / (y1 - y0)
    θz = z1 == z0 ? zero(T) : (z - z0) / (z1 - z0)

    c000 = vals[i, j, k]
    c100 = vals[i + 1, j, k]
    c010 = vals[i, j + 1, k]
    c110 = vals[i + 1, j + 1, k]
    c001 = vals[i, j, k + 1]
    c101 = vals[i + 1, j, k + 1]
    c011 = vals[i, j + 1, k + 1]
    c111 = vals[i + 1, j + 1, k + 1]

    c00 = _lerp(c000, c100, θx)
    c10 = _lerp(c010, c110, θx)
    c01 = _lerp(c001, c101, θx)
    c11 = _lerp(c011, c111, θx)
    c0 = _lerp(c00, c10, θy)
    c1 = _lerp(c01, c11, θy)
    return _lerp(c0, c1, θz)
end

function _interp_space(sf::SlabPhiField{1,T}, vals::AbstractVector{T}, x::T) where {T}
    return _interp1(sf.xyz[1], vals, x)
end

function _interp_space(sf::SlabPhiField{2,T}, vals::AbstractMatrix{T}, x::T, y::T) where {T}
    return _interp2(sf.xyz[1], sf.xyz[2], vals, x, y)
end

function _interp_space(sf::SlabPhiField{3,T}, vals::Array{T,3}, x::T, y::T, z::T) where {T}
    return _interp3(sf.xyz[1], sf.xyz[2], sf.xyz[3], vals, x, y, z)
end

function (sf::SlabPhiField{1,T})(x::T, t::T) where {T}
    θ = sf.t1 > sf.t0 ? clamp((t - sf.t0) / (sf.t1 - sf.t0), zero(T), one(T)) : zero(T)
    ϕn = _interp_space(sf, sf.phi_n, x)
    ϕp = _interp_space(sf, sf.phi_np1, x)
    return _lerp(ϕn, ϕp, θ)
end

function (sf::SlabPhiField{2,T})(x::T, y::T, t::T) where {T}
    θ = sf.t1 > sf.t0 ? clamp((t - sf.t0) / (sf.t1 - sf.t0), zero(T), one(T)) : zero(T)
    ϕn = _interp_space(sf, sf.phi_n, x, y)
    ϕp = _interp_space(sf, sf.phi_np1, x, y)
    return _lerp(ϕn, ϕp, θ)
end

function (sf::SlabPhiField{3,T})(x::T, y::T, z::T, t::T) where {T}
    θ = sf.t1 > sf.t0 ? clamp((t - sf.t0) / (sf.t1 - sf.t0), zero(T), one(T)) : zero(T)
    ϕn = _interp_space(sf, sf.phi_n, x, y, z)
    ϕp = _interp_space(sf, sf.phi_np1, x, y, z)
    return _lerp(ϕn, ϕp, θ)
end

function update_slab_field!(sf::SlabPhiField{N,T}, phi_n, phi_np1, t0::T, t1::T) where {N,T}
    size(phi_n) == size(sf.phi_n) || throw(DimensionMismatch("phi_n size mismatch"))
    size(phi_np1) == size(sf.phi_np1) || throw(DimensionMismatch("phi_np1 size mismatch"))
    sf.phi_n .= phi_n
    sf.phi_np1 .= phi_np1
    sf.t0 = t0
    sf.t1 = t1
    return sf
end
