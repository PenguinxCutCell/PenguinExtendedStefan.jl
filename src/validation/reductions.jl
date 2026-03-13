"""
    reduction_to_thermal_stefan(params)

Checks the concentration-decoupled thermal limit:
`m_liquidus = 0` and `k_partition = 1`.
"""
function reduction_to_thermal_stefan(params::ExtendedStefanParams; atol::Real=1e-12)
    ic = params.alloy_ic
    return abs(ic.m_liquidus) <= atol && abs(ic.k_partition - 1) <= atol
end

"""
    reduction_to_fixed_interface(params)

Heuristic check for near-fixed interface behavior in which `rhoL` is very large.
"""
function reduction_to_fixed_interface(params::ExtendedStefanParams; threshold::Real=1e8)
    return params.rhoL >= threshold
end
