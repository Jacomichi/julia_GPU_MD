using CUDA

include("md_struct.jl")

function ret_size(i::Int,para::Simu_Param)
    return ifelse(i <= para.N1, 1.0 ,para.σ_ratio)
end

function initial_lattice(pos,cell_num,simu_para::Simu_Param)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    (x,y) = (1,2)
    if i <= simu_para.N
        pos[x, i] = (simu_para.L - 1.0) * (floor((i - 1) / cell_num) / cell_num) + 1.#縦方向
        pos[y, i] = (simu_para.L - 1.0) * (((i - 1) % cell_num) / cell_num) + 1.#横方向
    end
    return nothing
end


#-------------------------------------------------------------------------------
#
#                       Initialize configuration
#
#
#-------------------------------------------------------------------------------



function initial_lattice(pos,L,cell_num,N)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    (x,y) = (1,2)
    if i <= N
        pos[x, i] = (L - 1.0) * (floor((i - 1) / cell_num) / cell_num) + 1.#縦方向
        pos[y, i] = (L - 1.0) * (((i - 1) % cell_num) / cell_num) + 1.#横方向
    end
    return nothing
end


function initial_random(pos,L,rnd)
    pos .= rand!(rnd) * L
end

function initialize_verlocity(velo)
    velo = CUDA.zeros(size(velo))
end


#-------------------------------------------------------------------------------
#
#                       list
#
#
#-------------------------------------------------------------------------------



function calc_grid_pos(pos_i::str2D,cell_size)
    grid_posx = ceil(Int64,pos_i.x/cell_size)
    grid_posy = ceil(Int64,pos_i.y/cell_size)
    grid_pos = str2D(grid_posx,grid_posy)
    return grid_pos
end


function calc_grid_address(grid_pos::str2D,Grid_num)
    return (grid_pos.y - 1) * Grid_num + grid_pos.x
end


function add_particle_to_cell(grid_pos,index,grid_counter,grid_cell,grid_para::Grid_Param)
    index = Int64(index)
    grid_address = calc_grid_address(grid_pos,grid_para.Grid_num)

    #CUDA.atomic_add!(ptr,val)の引数はptrと値(CUDA Cと同じ)なので、
    #CuDeviceArrayはpointer()を使って、pointerが取れることを使っている。
    #CuDeviceArrayはGPUから見た、CuArrayである。
    counter = CUDA.atomic_add!(pointer(grid_counter,grid_address), 1)

    counter = min(counter,grid_para.Max_particle_per_cell)

    grid_cell[(grid_address - 1) * grid_para.Max_particle_per_cell + counter + 1] = index
    return nothing
end


function update_grid(pos,grid_counter,grid_cell,grid_para::Grid_Param)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= size(pos)[2]
        pos_i = str2D(pos[1,index],pos[2,index])
        grid_pos = calc_grid_pos(pos_i,grid_para.cell_size)
        add_particle_to_cell(grid_pos,index,grid_counter,grid_cell,grid_para)
    end
    return nothing
end

function is_margin(velo,simu_para::Simu_Param,margin::Margin)
    v_max = maximum(hypot.(velo[1,:],velo[2,:]))
    margin.margin_length .-= 2.0 * v_max * simu_para.dt
    if sum(margin.margin_length) <= 0.0
        margin.margin_length .= margin.MARGIN
        return false
    end
        return true
end


#-------------------------------------------------------------------------------
#
#                                EOM
#
#
#-------------------------------------------------------------------------------



function calc_force_harmonic(pos1::str2D,pos2::str2D,p_size,simu_para::Simu_Param)
    force_x = zero(Float64)
    force_y = zero(Float64)

    distance = str2D(pos1.x - pos2.x,pos1.y - pos2.y)
    dist = fit_boundary(distance,simu_para.L)
    r = sqrt(dist.x * dist.x + dist.y * dist.y)
    if r < p_size
        force_x = dist.x/r * (1.0 - r/p_size)/p_size
        force_y = dist.y/r * (1.0 - r/p_size)/p_size
    end
    return str2D{Float64}(force_x,force_y)
end

function calc_force_WCA(pos1::str2D,pos2::str2D,p_size,simu_para::Simu_Param)
    force_x = zero(Float64)
    force_y = zero(Float64)

    distance = str2D(pos1.x - pos2.x,pos1.y - pos2.y)
    dist = fit_boundary(distance,simu_para.L)
    r = sqrt(dist.x * dist.x + dist.y * dist.y)
    if r < p_size * simu_para.WCA_cutoff
        invr = p_size/r
        invr² = invr * invr
        invR = invr² *  invr² *  invr²
        force_x = (24.0 * invR * (2.0 * invR - 1.0)) * invr² * dist.x
        force_y = (24.0 * invR * (2.0 * invR - 1.0)) * invr² * dist.y
    end
    return str2D{Float64}(force_x,force_y)
end


function calc_potential_harmonic(pos1::str2D,pos2::str2D,p_size,L)
    potential = 0.0

    distance = str2D(pos1.x - pos2.x,pos1.y - pos2.y)
    dist = fit_boundary(distance,L)
    r = sqrt(dist.x * dist.x + dist.y * dist.y)
    if r < p_size
        potential = (1.0 - r/p_size) * (1.0 - r/p_size)
    end
    return potential
end

function fit_boundary(dist::str2D,L)
    return str2D(dist.x - L * floor(dist.x/L + 0.5), dist.y - L * floor(dist.y/L + 0.5))
end

function fit_grid_cell(grid_pos,Grid_num)
    true_grid_posx = 0
    true_grid_posy = 0

    if grid_pos.x == 0
        true_grid_posx = Grid_num
    elseif grid_pos.x == Grid_num + 1
        true_grid_posx = 1
    else
        true_grid_posx = grid_pos.x
    end

    if grid_pos.y == 0
        true_grid_posy = Grid_num
    elseif grid_pos.y == Grid_num + 1
        true_grid_posy = 1
    else
        true_grid_posy = grid_pos.y
    end

    return str2D(true_grid_posx,true_grid_posy)
end

function collide_force_cell(grid_pos,index,s1,pos1::str2D,oldpos,grid_counter,grid_cells,simu_para::Simu_Param,grid_para::Grid_Param)
    #@cuprintf("%ld\n",index)
    new_grid_pos = fit_grid_cell(grid_pos,grid_para.Grid_num)
    grid_address = calc_grid_address(new_grid_pos,grid_para.Grid_num)
    particle_in_cell = grid_counter[grid_address]
    particle_in_cell = min(particle_in_cell,grid_para.Max_particle_per_cell)
    #@cuprintf("%ld\n",grid_address)
    force_x = zero(Float64)
    force_y = zero(Float64)
    for i = 1:particle_in_cell
        i_in_cell = (grid_address - 1) * grid_para.Max_particle_per_cell + i
        index2 = grid_cells[i_in_cell]

        if index2 != index

            pos2 = str2D(oldpos[1,index2],oldpos[2,index2])
            p_size = (s1 + ret_size(index2,simu_para))/2.0
            tmp_force = calc_force_WCA(pos1,pos2,p_size,simu_para)
            force_x += tmp_force.x
            force_y += tmp_force.y
        end
    end
    return str2D{Float64}(force_x,force_y)
end

function collide_potential_cell(grid_pos,index,s1,pos::str2D,oldpos,grid_counter,grid_cells,simu_para::Simu_Param,grid_para::Grid_Param)
    new_grid_pos = fit_grid_cell(grid_pos,grid_para.Grid_num)
    grid_address = calc_grid_address(new_grid_pos,grid_para.Grid_num)
    particle_in_cell = grid_counter[grid_address]
    particle_in_cell = min(particle_in_cell,grid_para.Max_particle_per_cell)

    potential = 0.0
    for i = 1:particle_in_cell
        index2 = grid_cell[(grid_address - 1) * grid_para.Max_particle_per_cell + i]

        if index2 != index
            pos2 = str2D(oldpos[1,index2],oldpos[2,index2])
            p_size = (s1 + ret_size(index2,simu_para))/2.0
            potential += calc_potential_harmonic(pos1,pos2,p_size,simu_para.L)
        end
    end
    return potential
end

function update_force!(pos,force,grid_counter,grid_cells,simu_para::Simu_Param,grid_para::Grid_Param)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= simu_para.N
        pos_i = str2D(pos[1,index],pos[2,index])
        force_x = zero(Float64)
        force_y = zero(Float64)
        grid_pos = calc_grid_pos(pos_i,grid_para.cell_size)
        s1 = ret_size(index,simu_para)
        for y in -1:1
            for x in -1:1
                search_grid_pos = str2D(grid_pos.x + x , grid_pos.y + y)
                tmp_force = collide_force_cell(search_grid_pos,index,s1,pos_i,pos,grid_counter,grid_cells,simu_para,grid_para)
                force_x += tmp_force.x
                force_y += tmp_force.y
            end
        end

        force[1,index] = force_x
        force[2,index] = force_y
    end

    return nothing
end


function update_velocity!(velo,force,rnd,simu_para::Simu_Param)
    velo .+= ( - simu_para.γ * simu_para.dt) .* velo .+ force .* simu_para.dt .+ simu_para.strength .* randn!(rnd)
end

function update_pos!(pos,velo,simu_para::Simu_Param)
    pos .+=  velo .* simu_para.dt
end

function check_boundary!(pos,simu_para::Simu_Param)
     pos .-= simu_para.L .* floor.(pos ./simu_para.L)
     return nothing
end

function calc_kinetic_energy(velo,simu_para::Simu_Param)
    energy = sum(velo .^ 2)/simu_para.N/2.0
    return energy
end


function calc_potential_energy(pos,energy,simu_para::Simu_Param,grid_para::Grid_Param,numthreads)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    tmp = CUDA.@cuDynamicSharedMem(Float64,numthreads)
    if index <= simu_para.N
        pos_i = str2D(pos[1,index],pos[2,index])
        potential = 0.0
        grid_pos = calc_grid_pos(pos_i,grid_para.cell_size)
        s1 = ret_size(index,simu_para)
        for y in -1:1
            for x in -1:1
                search_grid_pos = str2D(grid_pos.x + x.grid_pos.y + y)
                potential += collide_potential_cell(search_grid_pos,index,s1,pos_i,pos,grid_counter,grid_cells,para::Param)
            end
        end
        tmp[threadIdx().x] = potential
        CUDA.sync_threads()

        if threadIdx().x == 1
            sum = 0.0
            for i in 1:numthreads
                sum += temp[i]
            end
            CUDA.atomic_add!(pointer(energy), sum)
        end

    end
    return nothing
end
