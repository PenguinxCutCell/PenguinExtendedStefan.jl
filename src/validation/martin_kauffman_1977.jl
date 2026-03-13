struct MartinKauffmanReference{T}
    S::T
    tau::T
    Lambda::T
    A::T
    B::T
    alpha::T
    kappaT::T
end

@inline function _mk_erfc(x::T) where {T<:Real}
    z = abs(x)
    t = one(T) / (one(T) + T(0.5) * z)
    p = T(0.17087277)
    p = T(-0.82215223) + t * p
    p = T(1.48851587) + t * p
    p = T(-1.13520398) + t * p
    p = T(0.27886807) + t * p
    p = T(-0.18628806) + t * p
    p = T(0.09678418) + t * p
    p = T(0.37409196) + t * p
    p = T(1.00002368) + t * p
    τ = t * exp(-z * z - T(1.26551223) + t * p)
    return x >= zero(T) ? τ : (T(2) - τ)
end

function martin_kauffman_reference(; T::Type{<:Real}=Float64, kappaT::Real=1.0)
    TT = T
    return MartinKauffmanReference{TT}(
        TT(2.5),
        TT(0.1),
        TT(0.4),
        TT(0.90954),
        TT(0.44748),
        TT(0.19742),
        TT(kappaT),
    )
end

martin_kauffman_interface_position(t::T, pars::MartinKauffmanReference{T}) where {T} =
    t <= zero(T) ? zero(T) : 2 * pars.alpha * sqrt(pars.kappaT * t)

function martin_kauffman_fields(x::T, t::T, pars::MartinKauffmanReference{T}) where {T<:Real}
    t <= zero(T) && return (T=zero(T), C=zero(T), h=zero(T), eta=zero(T))

    η = x / (2 * sqrt(pars.kappaT * t))
    h = martin_kauffman_interface_position(t, pars)

    Tval = if x < h
        one(T) - pars.A * _mk_erfc(-η)
    else
        one(T) - pars.A * _mk_erfc(-pars.alpha)
    end

    Cval = if x < h
        one(T) - pars.B * _mk_erfc(-η / sqrt(pars.tau))
    else
        zero(T)
    end

    return (T=Tval, C=Cval, h=h, eta=η)
end

"""
Compatibility helper used by older tests/examples.
"""
function martin_kauffman_solution(
    x::T,
    t::T;
    κT::T,
    κS::T,
    A::T,
    B::T,
    α::T,
) where {T<:Real}
    pars = MartinKauffmanReference{T}(T(2.5), κS / κT, T(0.4), A, B, α, κT)
    return martin_kauffman_fields(x, t, pars)
end
