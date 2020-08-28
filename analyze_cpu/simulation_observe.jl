# ==============================================================================
#                        calculate variables
# ==============================================================================




function adjust_periodic(dr::Float64, sys_para::systemParameters)::Float64
    L_half = sys_para.L * 0.5
    dr = ifelse(dr > L_half,dr - sys_para.L,dr)
    dr = ifelse(dr < -L_half,dr + sys_para.L,dr)
    return dr
end



#calc g(r) from one configuration
#grのArrayを出力する
function calc_one_gr(r,dr,Brownians::Brownian,sys_para::systemParameters)
    count_max = ceil(Int,sys_para.L/2.0/dr)
    n = zeros(Int64,count_max)

    for i in 1:Brownians.N
        for j in 1:Brownians.N
            dx = Brownians.r[1, j] - Brownians.r[1, i]
            dy = Brownians.r[2, j] - Brownians.r[2, i]
            dx = adjust_periodic(dx,sys_para)
            dy = adjust_periodic(dy,sys_para)
            dist = hypot(dx,dy)
            if dist < r[count_max] && i != j
                R = ceil(Int,dist/dr)
                n[R] += 1
            end
        end
    end

    #ここってρ?ϕ?
    ρ = Brownians.N /(sys_para. L^2)
    para = 2.0 * π * ρ * Brownians.N * dr
    return n ./ r ./ para
    #=
    for i in 1:count_max
        para = 2.0 * π * ρ * r[i] * Brownians.N * dr
        gr[i] += n[i]/para
    end
    =#
end

function setting_gr(sys_para::systemParameters,dr)
    count_max = ceil(Int,sys_para.L/2.0/dr)
    r = [i*dr for i in 1:count_max]
    return r
end

#calc S(k)
function calc_sk(Brownians::Brownian,sys_para::systemParameters,k_min,k_max)
    k_count = ceil(Int,k_max/k_min)
    #println(typeof(k_count))
    Sk_xy = zeros(Float64,(k_count,k_count))

    for nx in 1:k_count
        kx = k_min*nx

        for ny in 1:k_count
            ky = k_min * ny

            a = 0.0
            b = 0.0
            for i in 1:Brownians.N
                a +=cos(Brownians.r[1,i] * kx + Brownians.r[2,i] * ky)
                b +=sin(Brownians.r[1,i] * kx + Brownians.r[2,i] * ky)
            end
            Sk_xy[nx,ny] += a * a + b * b
        end
    end

    Δk = k_min
    k_MAXcount = ceil(Int,sqrt(2.0) * k_max/Δk)
    sk = zeros(Float64,k_MAXcount)
    counter = zeros(Int64,k_MAXcount)
    for nx in 1:k_count
        for ny in 1:k_count
            distance = ceil(Int,hypot(nx,ny) * k_min / Δk)
            if distance < k_MAXcount
                sk[distance] += Sk_xy[nx,ny]/sys_para.N
                counter[distance] += 1
            end
        end
    end

    for i = 1:k_MAXcount
        if counter[i] != 0
            sk[i] = sk[i]/ counter[i]
        end
    end
    return sk
end

function calc_Fskt(IOSdata,k::Float64,Brownians::Brownian,sys_para::systemParameters)
    time = 0.0
    old_time = 0.0
    time_counter = 1
    sample_counter = 1
    time_arr = Vector{Float64}()
    Fskt = Vector{Float64}()
    Initial_config = zeros(Float64,2,Brownians.N)
    readlines_config_w_time(Brownians,time,Brownians.N,IOSdata)
    Initial_config = deepcopy(Brownians.r)
    while !eof(IOSdata)
        readlines_config_w_time(Brownians,time,Brownians.N,IOSdata)
        if old_time > time
            Initial_config = deepcopy(Brownians.r)
            time_counter = 1
            sample_counter += 1
        end
        for i in 1:Brownians.N
            dx = Brownians.r[1,i] - Initial_config[1,i]
            dy = Brownians.r[2,i] - Initial_config[2,i]
            dx = adjust_periodic(dx,sys_para)
            dy = adjust_periodic(dy,sys_para)
            a = cos(k * dx) + cos(k * dy)
            if counter ∉ keys(Fskt)
                push!(Fskt,a)
                push!(time_arr,time)
            else
                Fskt[counter] += a
            end
        end
        time_counter += 1
    end
    return time_arr,Fskt ./ sample_counter
end

function calc_velocity_correlation(r,dr,Brownians::Brownian,sys_para::systemParameters)
    count_max = ceil(sys_para.L/2.0/dr)
    xx,xy,yx,yy = 1,2,3,4
    velo_corr = zeros(Int64,4,count_max)

    for i in 1:Brownians.N
        for j in 1:Brownians.N
            dx = Brownians.r[1, j] - Brownians.r[1, i]
            dy = Brownians.r[2, j] - Brownians.r[2, i]
            dx = adjust_periodic(dx,sys_para)
            dy = adjust_periodic(dy,sys_para)
            r = hypot(dx,dy)
            if r < r[count_max] && i != j
                R = ceil(r/dr)
                velo_corr[xx,R] += Brownians.v[1, j] * Brownians.v[1, i]
                velo_corr[xy,R] += Brownians.v[1, j] * Brownians.v[2, i]
                velo_corr[yx,R] += Brownians.v[2, j] * Brownians.v[1, i]
                velo_corr[yy,R] += Brownians.v[2, j] * Brownians.v[2, i]
            end
        end
    end

    #ここってρ?ϕ?
    ρ = Brownians.N /(sys_para. L^2)
    para = 2.0 * π * ρ * Brownians.N * dr
    return velo_corr ./ r ./ para

end
