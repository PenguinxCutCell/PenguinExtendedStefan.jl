function _psi_vectors(Vn::AbstractVector{T}, Vn1::AbstractVector{T}, scheme) where {T}
    θ = PenguinDiffusion._theta_from_scheme(T, scheme)
    psip, psim = PenguinDiffusion._psi_functions(T, θ)
    ψp = Vector{T}(undef, length(Vn))
    ψm = Vector{T}(undef, length(Vn))
    @inbounds for i in eachindex(Vn)
        ψp[i] = convert(T, psip(Vn[i], Vn1[i]))
        ψm[i] = convert(T, psim(Vn[i], Vn1[i]))
    end
    return ψp, ψm
end

function _thermal_interface_fluxes(
    cache::ExtendedStefanCache{N,T},
    Tlω_new::AbstractVector{T},
    Tlγ_new::AbstractVector{T},
    Tsω_new::AbstractVector{T},
    Tsγ_new::AbstractVector{T},
    Tlω_prev::AbstractVector{T},
    Tlγ_prev::AbstractVector{T},
    Tsω_prev::AbstractVector{T},
    Tsγ_prev::AbstractVector{T},
    t::T,
    dt::T,
    scheme,
) where {N,T}
    model = cache.modelT
    cap_l = something(model.cap1_slab)
    ops_l = something(model.ops1_slab)
    cap_s = something(model.cap2_slab)
    ops_s = something(model.ops2_slab)

    θ = PenguinDiffusion._theta_from_scheme(T, scheme)
    _, _, Jl, Ll = PenguinDiffusion._weighted_core_ops(cap_l, ops_l, model.D1, t + θ * dt, model.coeff_mode)
    _, _, Js, Ls = PenguinDiffusion._weighted_core_ops(cap_s, ops_s, model.D2, t + θ * dt, model.coeff_mode)

    ψlp, ψlm = _psi_vectors(model.V1n, model.V1n1, scheme)
    ψsp, ψsm = _psi_vectors(model.V2n, model.V2n1, scheme)

    Tlω = ψlp .* Tlω_new .+ ψlm .* Tlω_prev
    Tlγ = ψlp .* Tlγ_new .+ ψlm .* Tlγ_prev
    Tsω = ψsp .* Tsω_new .+ ψsm .* Tsω_prev
    Tsγ = ψsp .* Tsγ_new .+ ψsm .* Tsγ_prev

    flux_l = Jl * Tlω + Ll * Tlγ
    flux_s = Js * Tsω + Ls * Tsγ
    return flux_l, flux_s
end

function stefan_speed!(
    v_nodes,
    frozen_mask::BitArray{N},
    cache::ExtendedStefanCache{N,T},
    Tlω_new::AbstractVector{T},
    Tlγ_new::AbstractVector{T},
    Tsω_new::AbstractVector{T},
    Tsγ_new::AbstractVector{T},
    Tlω_prev::AbstractVector{T},
    Tlγ_prev::AbstractVector{T},
    Tsω_prev::AbstractVector{T},
    Tsγ_prev::AbstractVector{T},
    rhoL,
    scheme,
    t::T,
    dt::T;
    tol::T=sqrt(eps(T)),
) where {N,T}
    cap_l = cache.modelT.cap1_slab
    cap_l === nothing && throw(ArgumentError("moving slab geometry is not initialized"))

    flux_l, flux_s = _thermal_interface_fluxes(
        cache,
        Tlω_new,
        Tlγ_new,
        Tsω_new,
        Tsγ_new,
        Tlω_prev,
        Tlγ_prev,
        Tsω_prev,
        Tsγ_prev,
        t,
        dt,
        scheme,
    )

    Γ = cap_l.buf.Γ
    vflat = vec(v_nodes)
    fflat = vec(frozen_mask)
    length(vflat) == length(Γ) || throw(DimensionMismatch("v_nodes must have $(length(Γ)) entries"))

    fill!(vflat, zero(T))
    fill!(fflat, false)

    rhoL_T = convert(T, rhoL)
    has_interface = false

    if N == 1
        sum_drive = zero(T)
        nactive = 0
        @inbounds for i in eachindex(vflat)
            γ = Γ[i]
            if isfinite(γ) && γ > tol
                sum_drive += flux_s[i] - flux_l[i]
                nactive += 1
                has_interface = true
            end
        end
        if nactive > 0
            # 1D moving-slab interface flux is integrated over the slab time
            # interval, so recover instantaneous speed with 1/dt scaling.
            v_if = (sum_drive / nactive) / (rhoL_T * dt)
            @inbounds for i in eachindex(vflat)
                γ = Γ[i]
                if isfinite(γ) && γ > tol
                    vflat[i] = v_if
                    fflat[i] = true
                end
            end
        end
    else
        @inbounds for i in eachindex(vflat)
            γ = Γ[i]
            if isfinite(γ) && γ > tol
                # For N>1, Γ carries the interface space-time measure from the moving
                # slab assembly, so dividing by Γ recovers the instantaneous speed.
                vflat[i] = (flux_s[i] - flux_l[i]) / (rhoL_T * γ)
                fflat[i] = true
                has_interface = true
            end
        end
    end

    if !has_interface
        @warn "No interface DOFs detected; Stefan speed set to zero for this step"
    end

    return v_nodes, frozen_mask
end

function thermal_stefan_residual(
    cache::ExtendedStefanCache{N,T},
    Tlω_new::AbstractVector{T},
    Tlγ_new::AbstractVector{T},
    Tsω_new::AbstractVector{T},
    Tsγ_new::AbstractVector{T},
    Tlω_prev::AbstractVector{T},
    Tlγ_prev::AbstractVector{T},
    Tsω_prev::AbstractVector{T},
    Tsγ_prev::AbstractVector{T},
    speed_full,
    rhoL,
    scheme,
    t::T,
    dt::T;
    tol::T=sqrt(eps(T)),
) where {N,T}
    cap_l = something(cache.modelT.cap1_slab)
    flux_l, flux_s = _thermal_interface_fluxes(
        cache,
        Tlω_new,
        Tlγ_new,
        Tsω_new,
        Tsγ_new,
        Tlω_prev,
        Tlγ_prev,
        Tsω_prev,
        Tsγ_prev,
        t,
        dt,
        scheme,
    )

    Γ = cap_l.buf.Γ
    v = vec(speed_full)
    r = zeros(T, length(Γ))
    maxr = zero(T)
    l2 = zero(T)
    n = 0
    rhoL_T = convert(T, rhoL)

    @inbounds for i in eachindex(Γ)
        γ = Γ[i]
        if !(isfinite(γ) && γ > tol)
            continue
        end
        # In 1D the integrated flux is converted with 1/dt, while in N>1 the
        # conversion uses Γ (space-time interface measure).
        drive = N == 1 ? (flux_s[i] - flux_l[i]) / dt : (flux_s[i] - flux_l[i]) / γ
        ri = rhoL_T * v[i] - drive
        r[i] = ri
        maxr = max(maxr, abs(ri))
        l2 += ri^2
        n += 1
    end

    l2 = n > 0 ? sqrt(l2 / n) : zero(T)
    return (maxabs=maxr, l2=l2, residual=r, count=n)
end

function solute_conservation_residual(
    cache::ExtendedStefanCache{N,T},
    Clω::AbstractVector{T},
    Clγ::AbstractVector{T},
    Csω::AbstractVector{T},
    Csγ::AbstractVector{T},
    speed_full,
    scheme,
    t::T,
    dt::T;
    tol::T=sqrt(eps(T)),
) where {N,T}
    model = cache.modelC
    cap_l = something(model.cap1_slab)
    ops_l = something(model.ops1_slab)
    cap_s = something(model.cap2_slab)
    ops_s = something(model.ops2_slab)

    θ = PenguinDiffusion._theta_from_scheme(T, scheme)
    _, _, Jl, Ll = PenguinDiffusion._weighted_core_ops(cap_l, ops_l, model.D1, t + θ * dt, model.coeff_mode)
    _, _, Js, Ls = PenguinDiffusion._weighted_core_ops(cap_s, ops_s, model.D2, t + θ * dt, model.coeff_mode)

    ql = Jl * Clω + Ll * Clγ
    qs = Js * Csω + Ls * Csγ

    Γ = cap_l.buf.Γ
    v = vec(speed_full)
    r = zeros(T, length(Γ))
    maxr = zero(T)
    l2 = zero(T)
    n = 0

    @inbounds for i in eachindex(Γ)
        γ = Γ[i]
        if !(isfinite(γ) && γ > tol)
            continue
        end
        ri = (qs[i] - ql[i]) - v[i] * (Clγ[i] - Csγ[i])
        r[i] = ri
        maxr = max(maxr, abs(ri))
        l2 += ri^2
        n += 1
    end

    l2 = n > 0 ? sqrt(l2 / n) : zero(T)
    return (maxabs=maxr, l2=l2, residual=r, count=n)
end
