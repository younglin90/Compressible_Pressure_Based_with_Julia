mutable struct controls
    Nx::UInt32
    Ny::UInt32
    Nz::UInt32

    Lx::Float64
    Ly::Float64
    Lz::Float64

    temporal_discretizationScheme::String
    spatial_discretizationScheme::String
    
    gravity::Array{Float64}

    p∞::Array{Float64}
	cᵥ::Array{Float64}
    γ::Array{Float64}
    b::Array{Float64}
	q::Array{Float64}
	

    left_p_BCtype
    left_u_BCtype
    left_v_BCtype
    left_T_BCtype
    left_Y_BCtype
    right_p_BCtype
    right_u_BCtype
    right_v_BCtype
    right_T_BCtype
    right_Y_BCtype
    bottom_p_BCtype
    bottom_u_BCtype
    bottom_v_BCtype
    bottom_T_BCtype
    bottom_Y_BCtype
    top_p_BCtype
    top_u_BCtype
    top_v_BCtype
    top_T_BCtype
    top_Y_BCtype
    left_p_BCValue
    left_u_BCValue
    left_v_BCValue
    left_T_BCValue
    left_Y_BCValue
    right_p_BCValue
    right_u_BCValue
    right_v_BCValue
    right_T_BCValue
    right_Y_BCValue
    bottom_p_BCValue
    bottom_u_BCValue
    bottom_v_BCValue
    bottom_T_BCValue
    bottom_Y_BCValue
    top_p_BCValue
    top_u_BCValue
    top_v_BCValue
    top_T_BCValue
    top_Y_BCValue

    
    realMaxIter::UInt32
    pseudoMaxIter::Vector{UInt32}
    pseudoMaxResidual::Float64

    corantNumber::Float64
    CFL::Float64
    
    Δt::Float64
    Lco::Float64
    Uco::Float64
    
    time::Float64
    realIter::UInt32
    pseudoIter::UInt32
    residual::Float64

    p::UInt32
    u::UInt32
    v::UInt32
    w::UInt32
    T::UInt32
    Y₁::UInt32
    ρ::UInt32
    Hₜ::UInt32
    c::UInt32
    Y₂::UInt32
    α₁::UInt32
    α₂::UInt32
    ∂ρ∂p::UInt32
    ∂ρ∂T::UInt32
    ∂ρ∂Y₁::UInt32
    ∂Hₜ∂p::UInt32
    ∂Hₜ∂T::UInt32
    ∂Hₜ∂Y₁::UInt32
    Δτ::UInt32
    Vᵣ::UInt32

    uⁿ::UInt32
    vⁿ::UInt32
    wⁿ::UInt32
    ρⁿ::UInt32
    pⁿ::UInt32
    Y₁ⁿ::UInt32
    Y₂ⁿ::UInt32
    α₁ⁿ::UInt32
    α₂ⁿ::UInt32
    Hₜⁿ::UInt32
    μ::UInt32
    

    uⁿ⁻¹::UInt32
    vⁿ⁻¹::UInt32
    wⁿ⁻¹::UInt32
    ρⁿ⁻¹::UInt32
    pⁿ⁻¹::UInt32
    Y₁ⁿ⁻¹::UInt32
    Y₂ⁿ⁻¹::UInt32
    α₁ⁿ⁻¹::UInt32
    α₂ⁿ⁻¹::UInt32
    Hₜⁿ⁻¹::UInt32
end
