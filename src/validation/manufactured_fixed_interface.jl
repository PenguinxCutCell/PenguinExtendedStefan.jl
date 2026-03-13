function manufactured_fixed_interface_1d_case(
    ;
    xΓ::Real=0.47,
    k_partition::Real=0.7,
    Tm::Real=0.0,
    m_liquidus::Real=-0.25,
    rho_cp_l::Real=1.0,
    rho_cp_s::Real=1.0,
    kappa_l::Real=1.0,
    kappa_s::Real=0.8,
    D_l::Real=0.12,
    D_s::Real=0.05,
)
    xg = float(xΓ)
    k = float(k_partition)
    Tm0 = float(Tm)
    mliq = float(m_liquidus)

    ClΓ(t) = 0.25 + 0.05 * sin(t)
    dClΓdt(t) = 0.05 * cos(t)
    CsΓ(t) = k * ClΓ(t)
    TΓ(t) = Tm0 + mliq * ClΓ(t)
    dTΓdt(t) = mliq * dClΓdt(t)

    Tl_exact(x, t) = TΓ(t) + 0.20 * (x - xg) + 0.08 * (x - xg)^2
    Ts_exact(x, t) = TΓ(t) - 0.15 * (x - xg) + 0.06 * (x - xg)^2
    Cl_exact(x, t) = ClΓ(t) + 0.10 * (x - xg) + 0.05 * (x - xg)^2
    Cs_exact(x, t) = CsΓ(t) - 0.08 * (x - xg) + 0.03 * (x - xg)^2

    source_T_l(x, t) = rho_cp_l * dTΓdt(t) - kappa_l * (2 * 0.08)
    source_T_s(x, t) = rho_cp_s * dTΓdt(t) - kappa_s * (2 * 0.06)
    source_C_l(x, t) = dClΓdt(t) - D_l * (2 * 0.05)
    source_C_s(x, t) = k * dClΓdt(t) - D_s * (2 * 0.03)

    return (
        xΓ=xg,
        k=k,
        Tm=Tm0,
        m_liquidus=mliq,
        phi0=(x) -> x - xg,
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

@inline function _fd_dt_2d(f, x, y, t; h=1e-6)
    return (f(x, y, t + h) - f(x, y, t - h)) / (2h)
end

@inline function _fd_lap_2d(f, x, y, t; h=5e-4)
    dxx = (f(x + h, y, t) - 2f(x, y, t) + f(x - h, y, t)) / (h^2)
    dyy = (f(x, y + h, t) - 2f(x, y, t) + f(x, y - h, t)) / (h^2)
    return dxx + dyy
end

function manufactured_fixed_interface_2d_case(
    ;
    center::Tuple{Real,Real}=(0.0, 0.0),
    radius::Real=0.35,
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
    R = float(radius)
    k = float(k_partition)
    Tm0 = float(Tm)
    mliq = float(m_liquidus)

    # Use a smooth polynomial interface function to avoid singular curvature terms.
    g(x, y) = (x - xc)^2 + (y - yc)^2 - R^2

    ClΓ(x, y, t) = 0.20 + 0.04 * sin(t) + 0.02 * (x - xc) - 0.01 * (y - yc)
    CsΓ(x, y, t) = k * ClΓ(x, y, t)
    TΓ(x, y, t) = Tm0 + mliq * ClΓ(x, y, t)

    aTl = 0.08
    aTs = -0.06
    aCl = 0.04
    # Match normal-flux balance at the interface for V=0: D_s*aCs = D_l*aCl.
    aCs = (D_l / D_s) * aCl

    Tl_exact(x, y, t) = TΓ(x, y, t) + aTl * g(x, y) + 0.02 * g(x, y)^2
    Ts_exact(x, y, t) = TΓ(x, y, t) + aTs * g(x, y) + 0.015 * g(x, y)^2
    Cl_exact(x, y, t) = ClΓ(x, y, t) + aCl * g(x, y) + 0.01 * g(x, y)^2
    Cs_exact(x, y, t) = CsΓ(x, y, t) + aCs * g(x, y) + 0.008 * g(x, y)^2

    source_T_l(x, y, t) = rho_cp_l * _fd_dt_2d(Tl_exact, x, y, t) - kappa_l * _fd_lap_2d(Tl_exact, x, y, t)
    source_T_s(x, y, t) = rho_cp_s * _fd_dt_2d(Ts_exact, x, y, t) - kappa_s * _fd_lap_2d(Ts_exact, x, y, t)
    source_C_l(x, y, t) = _fd_dt_2d(Cl_exact, x, y, t) - D_l * _fd_lap_2d(Cl_exact, x, y, t)
    source_C_s(x, y, t) = _fd_dt_2d(Cs_exact, x, y, t) - D_s * _fd_lap_2d(Cs_exact, x, y, t)

    return (
        center=(xc, yc),
        radius=R,
        k=k,
        Tm=Tm0,
        m_liquidus=mliq,
        phi0=(x, y) -> hypot(x - xc, y - yc) - R,
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
