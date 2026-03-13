struct ExtendedStefanParams{T,KT,DT,S1,S2,S3,S4,IC}
    rho_cp_l::T
    rho_cp_s::T
    kappa_l::KT
    kappa_s::KT
    D_l::DT
    D_s::DT
    rhoL::T
    source_T_l::S1
    source_T_s::S2
    source_C_l::S3
    source_C_s::S4
    alloy_ic::IC
end

function ExtendedStefanParams(
    rho_cp_l::T,
    rho_cp_s::T,
    rhoL::T,
    alloy_ic::AlloyEquilibrium;
    kappa_l=one(T),
    kappa_s=one(T),
    D_l=one(T),
    D_s=one(T),
    source_T_l=((args...) -> zero(T)),
    source_T_s=((args...) -> zero(T)),
    source_C_l=((args...) -> zero(T)),
    source_C_s=((args...) -> zero(T)),
) where {T<:Real}
    return ExtendedStefanParams{
        T,
        typeof(kappa_l),
        typeof(D_l),
        typeof(source_T_l),
        typeof(source_T_s),
        typeof(source_C_l),
        typeof(source_C_s),
        typeof(alloy_ic),
    }(
        rho_cp_l,
        rho_cp_s,
        kappa_l,
        kappa_s,
        D_l,
        D_s,
        rhoL,
        source_T_l,
        source_T_s,
        source_C_l,
        source_C_s,
        alloy_ic,
    )
end

struct ExtendedStefanOptions{T,LB,LI,ALG}
    scheme::Symbol
    ls_bc::LB
    ls_integrator::LI
    reinit::Bool
    reinit_every::Int
    extend_iters::Int
    extend_cfl::T
    interface_band::T
    coupling_max_iter::Int
    coupling_tol::T
    coupling_reltol::T
    coupling_damping::T
    alg::ALG
end

function ExtendedStefanOptions(
    ;
    scheme::Symbol=:BE,
    ls_bc=LevelSetMethods.NeumannBC(),
    ls_integrator=LevelSetMethods.RK3(),
    reinit::Bool=true,
    reinit_every::Int=1,
    extend_iters::Int=20,
    extend_cfl::Real=0.45,
    interface_band::Real=3.0,
    coupling_max_iter::Int=20,
    coupling_tol::Real=1e-10,
    coupling_reltol::Real=1e-10,
    coupling_damping::Real=0.5,
    alg=LinearSolve.KLUFactorization(),
)
    T = promote_type(
        typeof(float(extend_cfl)),
        typeof(float(interface_band)),
        typeof(float(coupling_tol)),
        typeof(float(coupling_reltol)),
        typeof(float(coupling_damping)),
    )
    if !(scheme === :BE || scheme === :CN)
        throw(ArgumentError("scheme must be :BE or :CN"))
    end
    return ExtendedStefanOptions{T,typeof(ls_bc),typeof(ls_integrator),typeof(alg)}(
        scheme,
        ls_bc,
        ls_integrator,
        reinit,
        reinit_every,
        extend_iters,
        convert(T, extend_cfl),
        convert(T, interface_band),
        coupling_max_iter,
        convert(T, coupling_tol),
        convert(T, coupling_reltol),
        convert(T, coupling_damping),
        alg,
    )
end

struct ExtendedStefanProblem{N,T,IR,OP,PA}
    grid::CartesianGrids.CartesianGrid{N,T}
    bcT::BorderConditions
    bcC::BorderConditions
    params::PA
    interface_rep::IR
    options::OP
end

function ExtendedStefanProblem(
    grid::CartesianGrids.CartesianGrid{N,T},
    bcT::BorderConditions,
    bcC::BorderConditions,
    params::ExtendedStefanParams,
    interface_rep::AbstractInterfaceRep{N,T},
    options::ExtendedStefanOptions=ExtendedStefanOptions(),
) where {N,T}
    return ExtendedStefanProblem{N,T,typeof(interface_rep),typeof(options),typeof(params)}(
        grid,
        bcT,
        bcC,
        params,
        interface_rep,
        options,
    )
end

mutable struct ExtendedStefanState{N,T,A}
    t::T
    Tlω::Vector{T}
    Tlγ::Vector{T}
    Tsω::Vector{T}
    Tsγ::Vector{T}
    Clω::Vector{T}
    Clγ::Vector{T}
    Csω::Vector{T}
    Csγ::Vector{T}
    speed_full::A
    frozen_mask::BitArray{N}
    logs::Dict{Symbol,Any}
end

mutable struct ExtendedStefanSolver{N,T,PR,ST,CA}
    problem::PR
    state::ST
    cache::CA
end
