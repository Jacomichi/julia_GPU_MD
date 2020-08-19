using CUDA
using Test, BenchmarkTools

include("md_host.jl")
include("md_struct.jl")

function main()
    N = 2^6
    ρ = 0.5
    T = 1.0
    γ = 1.0
    L = sqrt(N/ρ)
    σ_ratio = 1.0

    pos = CUDA.zeros(Float64,2,N)
    velo = CUDA.zeros(Float64,2,N)

    n_threads = ret_threadsize(N)
    n_blocks = ret_blocksize(N,n_threads)
    setting =Config(n_threads,n_blocks)

    initialize_Lattice(pos,L,N,setting)
    cell_num = floor(sqrt(N)) + 1
    print_lattice(Array(pos),cell_num)
    println("L=$(L)")

end


main()
