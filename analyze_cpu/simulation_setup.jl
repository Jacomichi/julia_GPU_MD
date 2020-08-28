# ==============================================================================
#                           Simulation parameters
# ==============================================================================
#outer constructorで書く
mutable struct systemParameters
    N::Int64
    ϕ::Float64
    L::Float64
    σ_ratio::Float64
end

function systemParameters(
    N::Int64,
    ρ::Float64,
    σ_ratio::Float64
)
    L = calc_L_rho(N,ρ)
    systemParameters(N, ρ, L,σ_ratio)
end

#binary用
function calc_L_binary(N::Int64,ϕ::Float64,σ_ratio::Float64)
    L = sqrt(N*π*(1.0 + σ_ratio^2)/(8.0*ϕ))
    return L
end

function calc_L_rho(N::Int64,ρ::Float64)
    L = sqrt(N/ρ)
    return L
end


# ==============================================================================
#                           Particle information
# ==============================================================================
mutable struct Brownian
    N::Int64
    r::Array{Float64,2}
    v::Array{Float64,2}
    σ::Array{Float64,1}
end

function Brownian(Num::Int64,σ_ratio::Float64)
    r = zeros(Float64, 2, Num)
    velocity = zeros(Float64, 2, Num) #ランダム分布
    half = Int64(Num/2)
    @assert half*2 == Num
    σ = vcat(fill(1.0,half),fill(σ_ratio,half))
    Brownian(Num, r, velocity, σ)
end

mutable struct Observable
    pressure::Array{Float64,1}
    potE::Array{Float64,1}
    timer::Array{Float64,1}
end

function Observable()
    pressure = zeros(Float64, 1)
    potE = zeros(Float64, 1)
    timer = zeros(Float64, 1)
    Observable(pressure, potE, timer)
end
