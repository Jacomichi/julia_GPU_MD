using CUDA
using Test, BenchmarkTools

include("md_host.jl")
include("md_struct.jl")
include("md_input.jl")
include("md_output.jl")


function main()
    id = return_ID(ARGS)
    #Simulation Parameter
    #コアを無駄にしないために、N=2^kの形にしておく
    dt = 1e-3

    ρ = 0.6
    N = 2^14
    N1 = Int64(2^10 / 2)
    σ_ratio = 1.0

    T = 1.0
    γ = 1.0
    #L = sqrt(N/ρ)

    simu_para = Simu_Param(dt,ρ,N,N1,σ_ratio,T,γ)

    #Grid parameter
    cell_size = 2.3
    grid_num = ceil(Int,simu_para.L / cell_size)
    println(grid_num)
    @assert grid_num > 3
    cell_size = simu_para.L / grid_num
    println(cell_size)
    max_particle_per_cell = 25

    grid_para = Grid_Param(cell_size,grid_num,max_particle_per_cell)

    MARGIN = 1.0
    margin = Margin(MARGIN)


    #simulation time parameter
    tau = Int64(1/simu_para.dt)
    equilbriation_tau = 100
    equlibriation_steps = equilbriation_tau * tau
    max_tau = 200
    max_steps = max_tau * tau

    #Simulation variable
    pos = CUDA.zeros(Float64,2,N)
    velo = CUDA.zeros(Float64,2,N)
    force = CUDA.zeros(Float64,2,N)
    rnd = CUDA.rand(Float64, 2,N)

    grid_counter = CUDA.zeros(Int64,grid_para.Grid_num * grid_para.Grid_num)
    grid_cells =CUDA.zeros(Int64,grid_para.Max_particle_per_cell * grid_para.Grid_num * grid_para.Grid_num)

    #GPU setting
    n_threads = ret_threadsize(N)
    n_blocks = ret_blocksize(N,n_threads)
    @assert N == n_threads * n_blocks
    setting =　Config(n_threads,n_blocks)
    println(n_threads)
    println(n_blocks)

    #make_outputfile
    #file_dir = create_file_harmonic_config(simu_para,id)
    file_dir = create_file_config("WCA_config","WCA",simu_para,id)
    IOwritefile = open_writefile(file_dir)

    #Initialize system
    #initialize_System_random(pos,velo,rnd,simu_para)
    initialize_System_Lattice(pos,velo,simu_para,setting)

    #Equilibration
    println("start equilibration")
    for t in 1:equilbriation_tau

        for s in 1:tau
            check_grid(pos,velo,grid_counter,grid_cells,simu_para,grid_para,margin,setting)
            collide(pos,velo,force,rnd,grid_counter,grid_cells,simu_para,grid_para,setting)
        end
    end

    println("finish equlibriation")


    #main simulation
    println("start main simulation")
    for t in 1:max_tau

        for s in 1:tau
            check_grid(pos,velo,grid_counter,grid_cells,simu_para,grid_para,margin,setting)
            collide(pos,velo,force,rnd,grid_counter,grid_cells,simu_para,grid_para,setting)
        end

        host_pos = Array(pos)
        output_position(pos, IOwritefile)
    end

    println("finish main simulation")

    kinE = calc_kinetic_energy(velo,simu_para)

    println("final kinetic energy is ",kinE )


end


main()
