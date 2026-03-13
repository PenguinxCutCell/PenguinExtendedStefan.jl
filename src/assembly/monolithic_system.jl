function assemble_monolithic_system!(
    cache::ExtendedStefanCache{N,T},
    phi_n,
    phi_np1,
    Vlag_full,
    Tlω_prev::AbstractVector{T},
    Tlγ_prev::AbstractVector{T},
    Tsω_prev::AbstractVector{T},
    Tsγ_prev::AbstractVector{T},
    Clω_prev::AbstractVector{T},
    Clγ_prev::AbstractVector{T},
    Csω_prev::AbstractVector{T},
    Csγ_prev::AbstractVector{T},
    t::T,
    dt::T,
    prob::ExtendedStefanProblem{N,T},
) where {N,T}
    dt > zero(T) || throw(ArgumentError("dt must be positive"))

    update_slab_field!(cache.slab, phi_n, phi_np1, t, t + dt)
    PenguinDiffusion._build_moving_slab!(cache.modelT, t, dt)
    PenguinDiffusion._build_moving_slab!(cache.modelC, t, dt)

    capTl = something(cache.modelT.cap1_slab)
    opsTl = something(cache.modelT.ops1_slab)
    capTs = something(cache.modelT.cap2_slab)
    opsTs = something(cache.modelT.ops2_slab)

    capCl = something(cache.modelC.cap1_slab)
    opsCl = something(cache.modelC.ops1_slab)
    capCs = something(cache.modelC.cap2_slab)
    opsCs = something(cache.modelC.ops2_slab)

    nt = cache.layout.nt
    lay = cache.layout
    nsys = 8 * nt

    θ = PenguinDiffusion._theta_from_scheme(T, prob.options.scheme)
    psip, psim = PenguinDiffusion._psi_functions(T, θ)

    KTl, CTl, _, _ = PenguinDiffusion._weighted_core_ops(capTl, opsTl, cache.modelT.D1, t + θ * dt, cache.modelT.coeff_mode)
    KTs, CTs, _, _ = PenguinDiffusion._weighted_core_ops(capTs, opsTs, cache.modelT.D2, t + θ * dt, cache.modelT.coeff_mode)
    KCl, CCl, JCl, LCl = PenguinDiffusion._weighted_core_ops(capCl, opsCl, cache.modelC.D1, t + θ * dt, cache.modelC.coeff_mode)
    KCs, CCs, JCs, LCs = PenguinDiffusion._weighted_core_ops(capCs, opsCs, cache.modelC.D2, t + θ * dt, cache.modelC.coeff_mode)

    fTl_n, fTs_n = PenguinDiffusion._source_values_diph(capTl, cache.modelT.source1, capTs, cache.modelT.source2, t)
    fTl_n1, fTs_n1 = PenguinDiffusion._source_values_diph(capTl, cache.modelT.source1, capTs, cache.modelT.source2, t + dt)
    fCl_n, fCs_n = PenguinDiffusion._source_values_diph(capCl, cache.modelC.source1, capCs, cache.modelC.source2, t)
    fCl_n1, fCs_n1 = PenguinDiffusion._source_values_diph(capCl, cache.modelC.source1, capCs, cache.modelC.source2, t + dt)

    MTLn = spdiagm(0 => cache.modelT.V1n)
    MTLn1 = spdiagm(0 => cache.modelT.V1n1)
    MTSn = spdiagm(0 => cache.modelT.V2n)
    MTSn1 = spdiagm(0 => cache.modelT.V2n1)

    MCLn = spdiagm(0 => cache.modelC.V1n)
    MCLn1 = spdiagm(0 => cache.modelC.V1n1)
    MCSn = spdiagm(0 => cache.modelC.V2n)
    MCSn1 = spdiagm(0 => cache.modelC.V2n1)

    ψTLp = T[psip(cache.modelT.V1n[i], cache.modelT.V1n1[i]) for i in 1:nt]
    ψTLm = T[psim(cache.modelT.V1n[i], cache.modelT.V1n1[i]) for i in 1:nt]
    ψTSp = T[psip(cache.modelT.V2n[i], cache.modelT.V2n1[i]) for i in 1:nt]
    ψTSm = T[psim(cache.modelT.V2n[i], cache.modelT.V2n1[i]) for i in 1:nt]

    ψCLp = T[psip(cache.modelC.V1n[i], cache.modelC.V1n1[i]) for i in 1:nt]
    ψCLm = T[psim(cache.modelC.V1n[i], cache.modelC.V1n1[i]) for i in 1:nt]
    ψCSp = T[psip(cache.modelC.V2n[i], cache.modelC.V2n1[i]) for i in 1:nt]
    ψCSm = T[psim(cache.modelC.V2n[i], cache.modelC.V2n1[i]) for i in 1:nt]

    ΨTLp = spdiagm(0 => ψTLp)
    ΨTLm = spdiagm(0 => ψTLm)
    ΨTSp = spdiagm(0 => ψTSp)
    ΨTSm = spdiagm(0 => ψTSm)

    ΨCLp = spdiagm(0 => ψCLp)
    ΨCLm = spdiagm(0 => ψCLm)
    ΨCSp = spdiagm(0 => ψCSp)
    ΨCSm = spdiagm(0 => ψCSm)

    ρcp_l = convert(T, prob.params.rho_cp_l)
    ρcp_s = convert(T, prob.params.rho_cp_s)

    A_Tlω_Tlω = ρcp_l * MTLn1 + θ * (KTl * ΨTLp)
    A_Tlω_Tlγ = -ρcp_l * (MTLn1 - MTLn) + θ * (CTl * ΨTLp)

    A_Tsω_Tsω = ρcp_s * MTSn1 + θ * (KTs * ΨTSp)
    A_Tsω_Tsγ = -ρcp_s * (MTSn1 - MTSn) + θ * (CTs * ΨTSp)

    A_Clω_Clω = MCLn1 + θ * (KCl * ΨCLp)
    A_Clω_Clγ = -(MCLn1 - MCLn) + θ * (CCl * ΨCLp)

    A_Csω_Csω = MCSn1 + θ * (KCs * ΨCSp)
    A_Csω_Csγ = -(MCSn1 - MCSn) + θ * (CCs * ΨCSp)

    b_Tlω = (ρcp_l * MTLn - (one(T) - θ) * (KTl * ΨTLm)) * Tlω_prev
    b_Tlω .-= (one(T) - θ) .* ((CTl * ΨTLm) * Tlγ_prev)
    b_Tlω .+= θ .* (capTl.V * fTl_n1) .+ (one(T) - θ) .* (capTl.V * fTl_n)

    b_Tsω = (ρcp_s * MTSn - (one(T) - θ) * (KTs * ΨTSm)) * Tsω_prev
    b_Tsω .-= (one(T) - θ) .* ((CTs * ΨTSm) * Tsγ_prev)
    b_Tsω .+= θ .* (capTs.V * fTs_n1) .+ (one(T) - θ) .* (capTs.V * fTs_n)

    b_Clω = (MCLn - (one(T) - θ) * (KCl * ΨCLm)) * Clω_prev
    b_Clω .-= (one(T) - θ) .* ((CCl * ΨCLm) * Clγ_prev)
    b_Clω .+= θ .* (capCl.V * fCl_n1) .+ (one(T) - θ) .* (capCl.V * fCl_n)

    b_Csω = (MCSn - (one(T) - θ) * (KCs * ΨCSm)) * Csω_prev
    b_Csω .-= (one(T) - θ) .* ((CCs * ΨCSm) * Csγ_prev)
    b_Csω .+= θ .* (capCs.V * fCs_n1) .+ (one(T) - θ) .* (capCs.V * fCs_n)

    A = spzeros(T, nsys, nsys)
    b = zeros(T, nsys)

    PenguinDiffusion._insert_block!(A, lay.Tlω, lay.Tlω, A_Tlω_Tlω)
    PenguinDiffusion._insert_block!(A, lay.Tlω, lay.Tlγ, A_Tlω_Tlγ)
    PenguinDiffusion._insert_block!(A, lay.Tsω, lay.Tsω, A_Tsω_Tsω)
    PenguinDiffusion._insert_block!(A, lay.Tsω, lay.Tsγ, A_Tsω_Tsγ)
    PenguinDiffusion._insert_block!(A, lay.Clω, lay.Clω, A_Clω_Clω)
    PenguinDiffusion._insert_block!(A, lay.Clω, lay.Clγ, A_Clω_Clγ)
    PenguinDiffusion._insert_block!(A, lay.Csω, lay.Csω, A_Csω_Csω)
    PenguinDiffusion._insert_block!(A, lay.Csω, lay.Csγ, A_Csω_Csγ)

    PenguinDiffusion._insert_vec!(b, lay.Tlω, b_Tlω)
    PenguinDiffusion._insert_vec!(b, lay.Tsω, b_Tsω)
    PenguinDiffusion._insert_vec!(b, lay.Clω, b_Clω)
    PenguinDiffusion._insert_vec!(b, lay.Csω, b_Csω)

    apply_interface_rows!(
        A,
        b,
        cache,
        prob,
        capCl,
        capCs,
        JCl,
        LCl,
        JCs,
        LCs,
        Vlag_full,
        t + dt,
    )

    cache.sys.A = A
    cache.sys.b = b
    length(cache.sys.x) == nsys || (cache.sys.x = zeros(T, nsys))
    cache.sys.cache = nothing

    layTl = UnknownLayout(nt, (ω=lay.Tlω,))
    layTs = UnknownLayout(nt, (ω=lay.Tsω,))
    layCl = UnknownLayout(nt, (ω=lay.Clω,))
    layCs = UnknownLayout(nt, (ω=lay.Csω,))

    apply_box_bc_mono!(cache.sys.A, cache.sys.b, capTl, opsTl, cache.modelT.D1, cache.modelT.bc_border; t=t + θ * dt, layout=layTl)
    apply_box_bc_mono!(cache.sys.A, cache.sys.b, capTs, opsTs, cache.modelT.D2, cache.modelT.bc_border; t=t + θ * dt, layout=layTs)
    apply_box_bc_mono!(cache.sys.A, cache.sys.b, capCl, opsCl, cache.modelC.D1, cache.modelC.bc_border; t=t + θ * dt, layout=layCl)
    apply_box_bc_mono!(cache.sys.A, cache.sys.b, capCs, opsCs, cache.modelC.D2, cache.modelC.bc_border; t=t + θ * dt, layout=layCs)

    active_rows = _active_rows_extended(capCl, capCs, lay)
    cache.sys.A, cache.sys.b = PenguinDiffusion._apply_row_identity_constraints!(cache.sys.A, cache.sys.b, active_rows)

    return cache.sys
end

function solve_monolithic_system!(
    cache::ExtendedStefanCache{N,T},
    phi_n,
    phi_np1,
    Vlag_full,
    Tlω_prev::AbstractVector{T},
    Tlγ_prev::AbstractVector{T},
    Tsω_prev::AbstractVector{T},
    Tsγ_prev::AbstractVector{T},
    Clω_prev::AbstractVector{T},
    Clγ_prev::AbstractVector{T},
    Csω_prev::AbstractVector{T},
    Csγ_prev::AbstractVector{T},
    t::T,
    dt::T,
    prob::ExtendedStefanProblem{N,T};
    method::Symbol=:linearsolve,
    kwargs...,
) where {N,T}
    sys = assemble_monolithic_system!(
        cache,
        phi_n,
        phi_np1,
        Vlag_full,
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
        prob,
    )

    if method === :linearsolve
        solve!(sys; method=:linearsolve, reuse_factorization=false, alg=prob.options.alg, kwargs...)
    else
        solve!(sys; method=method, reuse_factorization=false, kwargs...)
    end

    lay = cache.layout
    Tlω_new = Vector{T}(sys.x[lay.Tlω])
    Tlγ_new = Vector{T}(sys.x[lay.Tlγ])
    Tsω_new = Vector{T}(sys.x[lay.Tsω])
    Tsγ_new = Vector{T}(sys.x[lay.Tsγ])
    Clω_new = Vector{T}(sys.x[lay.Clω])
    Clγ_new = Vector{T}(sys.x[lay.Clγ])
    Csω_new = Vector{T}(sys.x[lay.Csω])
    Csγ_new = Vector{T}(sys.x[lay.Csγ])

    uprev = _as_prev_full_extended(
        cache,
        Tlω_prev,
        Tlγ_prev,
        Tsω_prev,
        Tsγ_prev,
        Clω_prev,
        Clγ_prev,
        Csω_prev,
        Csγ_prev,
    )

    return Tlω_new, Tlγ_new, Tsω_new, Tsγ_new, Clω_new, Clγ_new, Csω_new, Csγ_new, uprev, sys
end
