using CUDA
using Test, BenchmarkTools
using Random: rand!,randn!


function ret_blocksize(N::Int64,numthreads::Int64)
    return ceil.(Int,N/numthreads)
end

function ret_threadsize(N::Int64)
    if N <= 1024
        return 2^(floor(Int,log(2,N)))
    else
        return 1024
    end
end

function update_pos(pos::CuArray{Float64,2},velo::CuArray{Float64,2},dt::Float64)
    pos = velo * dt
end

function Langevin(
    pos::CuArray{Float64,2},velo::CuArray{Float64,2},rnd::CuArray{Float64,2},
    dt::Float64,γ::Float64,T::Float64
)
    noise_strength = 2.0 * γ * T / dt
    velo .+= (γ .* velo .+ sqrt(noise_strength) * randn!(rnd)) .* simu_para.dt
end

function check_BC(pos::Float64, L::Float64)::Float64
    pos = ifelse(pos < 0.0,L, pos)
    pos = ifelse(pos > L,pos - L,pos)
    return pos
end

function initial_lattice(pos,L,cell_num,N)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    (x,y) = (1,2)
    if i <= N
        pos[x, i] = (L - 6.0) * (floor((i - 1) / cell_num) / cell_num) + 3#縦方向
        pos[y, i] = (L - 6.0) * (((i - 1) % cell_num) / cell_num) + 3#横方向
    end
    return nothing
end

function print_lattice(arr,cellnum)
    count = 0
    for j in 1:size(arr)[2]
        print("(")
        for i in 1:size(arr)[1]
            print("$(round(arr[i,j],digits=2)) ")
        end
        print(")")
        count += 1
        if count == cellnum
            count = 0
            println("")
        end
    end
    println("")
end

function initialize_Lattice(pos,L,N,setting)
    cell_num = floor(sqrt(N)) + 1
    @cuda threads=setting.numthreads blocks=setting.numblocks initial_lattice(pos,L,cell_num,N)
end
