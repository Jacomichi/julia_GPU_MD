using CUDA
using Test, BenchmarkTools
using Random: rand!,randn!

include("md_kernel.jl")
include("md_struct.jl")

#-------------------------------------------------------------------------------
#
#                       Setting GPU
#
#
#-------------------------------------------------------------------------------



function ret_blocksize(N::Int64,numthreads::Int64)
    return ceil.(Int,N/numthreads)
end

#ここでthreadのサイズを2^8をmaxにしているのは、
#今のコードだと、これ以上はレジスタを食い尽くして、止まってしまうため。
function ret_threadsize(N::Int64)
    if N <= 2^8
        return 2^(floor(Int,log(2,N)))
    else
        return 2^8
    end
end



#-------------------------------------------------------------------------------
#
#                       Initialization
#
#-------------------------------------------------------------------------------


function initialize_System_Lattice(pos,velo,simu_para::Simu_Param,setting::Config)
    cell_num = floor(sqrt(simu_para.N)) + 1
    @cuda threads=setting.numthreads blocks=setting.numblocks initial_lattice(pos,cell_num,simu_para)
    initialize_verlocity(velo)
end

function initialize_System_random(pos,velo,rnd,simu_para::Simu_Param)
    initial_random(pos,simu_para.L,rnd)
    initialize_verlocity(velo)
end


#-------------------------------------------------------------------------------
#
#                       Simulation
#
#-------------------------------------------------------------------------------

function check_grid(pos,velo,grid_counter,grid_cell,simu_para::Simu_Param,grid_para::Grid_Param,margin::Margin,setting::Config)
    if !is_margin(velo,simu_para,margin)
        update_grid(pos,grid_counter,grid_cell,grid_para,setting)
    end
    return nothing
end

function update_grid(pos,grid_counter,grid_cell,grid_para::Grid_Param,setting::Config)
    grid_counter .= CUDA.zeros(Int64,length(grid_counter))
    @cuda threads=setting.numthreads blocks=setting.numblocks update_grid(pos,grid_counter,grid_cell,grid_para)
end

function collide(pos,velo,force,rnd,grid_counter,grid_cells,simu_para::Simu_Param,grid_para::Grid_Param,setting::Config)
    update_pos!(pos,velo,simu_para)
    check_boundary!(pos,simu_para)
    @cuda threads=setting.numthreads blocks=setting.numblocks update_force!(pos,force,grid_counter,grid_cells,simu_para,grid_para)
    update_velocity!(velo,force,rnd,simu_para)
end

function calc_potential(pos,energy,simu_para::Simu_Param,grid_para::Grid_Param,setting::Config)
    energy = CUDA.zeros(Float64,1)
    @cuda threads=setting.numthreads blocks=setting.numblocks calc_potential_energy(pos,energy,simu_para,grid_para,numthreads,setting.numthreads)
    return energy[1]
end


#function run_simulation(pos,velo,force,rnd,grid_counter,grid_cells,simu_para::Simu_Param,grid_para::Grid_Param,setting::Config)

#-------------------------------------------------------------------------------
#
#                       Unitlity tools
#
#-------------------------------------------------------------------------------


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
