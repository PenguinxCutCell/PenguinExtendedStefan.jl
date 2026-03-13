struct ExtendedLayout
    nt::Int
    Tlω::UnitRange{Int}
    Tlγ::UnitRange{Int}
    Tsω::UnitRange{Int}
    Tsγ::UnitRange{Int}
    Clω::UnitRange{Int}
    Clγ::UnitRange{Int}
    Csω::UnitRange{Int}
    Csγ::UnitRange{Int}
end

function ExtendedLayout(nt::Int)
    return ExtendedLayout(
        nt,
        1:nt,
        (nt + 1):(2 * nt),
        (2 * nt + 1):(3 * nt),
        (3 * nt + 1):(4 * nt),
        (4 * nt + 1):(5 * nt),
        (5 * nt + 1):(6 * nt),
        (6 * nt + 1):(7 * nt),
        (7 * nt + 1):(8 * nt),
    )
end

mutable struct ExtendedStefanCache{N,T,SF,MT,MC,SY}
    slab::SF
    modelT::MT
    modelC::MC
    sys::SY
    layout::ExtendedLayout
    metrics::Dict{Symbol,Any}
end

function build_cache(prob::ExtendedStefanProblem{N,T}) where {N,T}
    ϕ = phi_values(prob.interface_rep)
    xyz = ntuple(d -> collect(T, CartesianGrids.grid1d(prob.grid, d)), N)
    slab = SlabPhiField{N,T,typeof(copy(ϕ))}(xyz, copy(ϕ), copy(ϕ), zero(T), one(T))

    modelT = PenguinDiffusion.MovingDiffusionModelDiph(
        prob.grid,
        slab,
        prob.params.kappa_l,
        prob.params.kappa_s;
        source=(prob.params.source_T_l, prob.params.source_T_s),
        bc_border=prob.bcT,
        ic=nothing,
        coeff_mode=:harmonic,
        geom_method=:vofijul,
    )

    modelC = PenguinDiffusion.MovingDiffusionModelDiph(
        prob.grid,
        slab,
        prob.params.D_l,
        prob.params.D_s;
        source=(prob.params.source_C_l, prob.params.source_C_s),
        bc_border=prob.bcC,
        ic=nothing,
        coeff_mode=:harmonic,
        geom_method=:vofijul,
    )

    nt = prod(prob.grid.n)
    lay = ExtendedLayout(nt)
    nsys = 8 * nt
    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys))

    metrics = Dict{Symbol,Any}(
        :solute_residual_norm => T(Inf),
        :stefan_residual_norm => T(Inf),
        :speed_increment_norm => T(Inf),
        :max_speed => zero(T),
        :nonlinear_iters => 0,
    )

    return ExtendedStefanCache{N,T,typeof(slab),typeof(modelT),typeof(modelC),typeof(sys)}(
        slab,
        modelT,
        modelC,
        sys,
        lay,
        metrics,
    )
end

function _as_prev_full_extended(
    cache::ExtendedStefanCache{N,T},
    Tlω_prev::AbstractVector{T},
    Tlγ_prev::AbstractVector{T},
    Tsω_prev::AbstractVector{T},
    Tsγ_prev::AbstractVector{T},
    Clω_prev::AbstractVector{T},
    Clγ_prev::AbstractVector{T},
    Csω_prev::AbstractVector{T},
    Csγ_prev::AbstractVector{T},
) where {N,T}
    nt = cache.layout.nt
    length(Tlω_prev) == nt || throw(DimensionMismatch("Tlω_prev length must be $nt"))
    length(Tlγ_prev) == nt || throw(DimensionMismatch("Tlγ_prev length must be $nt"))
    length(Tsω_prev) == nt || throw(DimensionMismatch("Tsω_prev length must be $nt"))
    length(Tsγ_prev) == nt || throw(DimensionMismatch("Tsγ_prev length must be $nt"))
    length(Clω_prev) == nt || throw(DimensionMismatch("Clω_prev length must be $nt"))
    length(Clγ_prev) == nt || throw(DimensionMismatch("Clγ_prev length must be $nt"))
    length(Csω_prev) == nt || throw(DimensionMismatch("Csω_prev length must be $nt"))
    length(Csγ_prev) == nt || throw(DimensionMismatch("Csγ_prev length must be $nt"))

    lay = cache.layout
    u = zeros(T, 8 * nt)
    u[lay.Tlω] .= Tlω_prev
    u[lay.Tlγ] .= Tlγ_prev
    u[lay.Tsω] .= Tsω_prev
    u[lay.Tsγ] .= Tsγ_prev
    u[lay.Clω] .= Clω_prev
    u[lay.Clγ] .= Clγ_prev
    u[lay.Csω] .= Csω_prev
    u[lay.Csγ] .= Csγ_prev
    return u
end
