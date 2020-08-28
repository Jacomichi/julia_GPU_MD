using CUDA

struct str2D{T<:Real}
    x::T
    y::T
end

struct Simu_Param
    dt::Float64
    ρ::Float64
    N::Int64
    N1::Int64
    L::Float64
    σ_ratio::Float64
    strength::Float64
    T::Float64
    γ::Float64
    WCA_cutoff::Float64
end

function Simu_Param(dt,ρ,N,N1,σ_ratio,T,γ)
    L = sqrt(N/ρ)
    strength = sqrt(2.0 * γ * T * dt)
    WCA_cutoff = 1.122462048309373
    Simu_Param(dt,ρ,N,N1,L,σ_ratio,strength,T,γ,WCA_cutoff)
end

struct Grid_Param
    cell_size::Float64
    Grid_num::Int64
    Max_particle_per_cell::Int64
end

struct Margin
    MARGIN::Float64
    margin_length::CuArray{Float64,1}
end

function Margin(MARGIN::Float64)
    Margin(MARGIN,CUDA.zeros(Float64,1))
end

struct Config
    numthreads::Int64
    numblocks::Int64
end
