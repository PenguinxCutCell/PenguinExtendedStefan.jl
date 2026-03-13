"""
    manufactured_alloy_planar_1d(x, t; s0, V, Tm, m, k, Cl_if)

Simple manufactured planar moving-interface state used for smoke/consistency tests.
"""
function manufactured_alloy_planar_1d(
    x::T,
    t::T;
    s0::T,
    V::T,
    Tm::T,
    m::T,
    k::T,
    Cl_if::T,
) where {T<:Real}
    s = s0 + V * t
    Tl = x <= s ? (Tm + m * Cl_if + (s - x)) : (Tm + m * Cl_if)
    Ts = x >= s ? (Tm + m * Cl_if + (x - s)) : (Tm + m * Cl_if)
    Cl = x <= s ? Cl_if : zero(T)
    Cs = x >= s ? k * Cl_if : zero(T)
    return (Tl=Tl, Ts=Ts, Cl=Cl, Cs=Cs, s=s)
end
