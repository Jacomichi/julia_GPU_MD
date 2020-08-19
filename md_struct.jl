using CUDA


struct Config
    numthreads::Int64
    numblocks::Int64
end

struct systemParameters
    N::Int64
    ρ::Float64
    L::Float64
    T::Float64
    γ::Float64
    σ_ratio::Float64
end

function systemParameters(
    N::Int64,
    ρ::Float64,
    T::Float64,
    γ::Float64,
    σ_ratio::Float64
)
    L = sqrt(N/ρ)
    systemParameters(N, ρ, L, T, γ,σ_ratio)
end


struct simulationParameters
    dt::Float64
    Tmax::Int64
    Tsteps::Int64
    outFreq::Float64
end

function simulationParameters(dt::Float64, Tmax::Int64, outFreq::Float64)
    Tsteps = floor(Tmax / dt)
    simulationParameters(dt, Tmax, Tsteps, outFreq)
end


mutable struct Brownian
    N::Int64
    r::CuArray{Float64,2}
    v::CuArray{Float64,2}
    rnd::CuArray{Float64,2}
end

function Brownian(Num::Int64,L::Float64)
    r = L * CUDA.rand(Float64, 2, Num) #ランダム分布
    v = CUDA.zeros(Float64, size(r)...) #速度は全て0
    rnd = CuArray{Float64}(undef, size(r)...)
    Brownian(Num, r, v,rnd)
end
