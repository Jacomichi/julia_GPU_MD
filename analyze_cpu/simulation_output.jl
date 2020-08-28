# ==============================================================================
#                           output data
# ==============================================================================


include("simulation_setup.jl")


function open_writefile(openfile::String)::IOStream
    outfile = open(openfile,"w")
end

function output_position(Brownians::Brownian, outfile::IOStream)
    for i = 1:Brownians.N
        println(outfile,"$(Brownians.r[1,i]) $(Brownians.r[2,i])")
    end
end


function output_config_data(Brownians::Brownian, outfile::IOStream)
    for i = 1:Brownians.N
        println(outfile,"$(Brownians.r[1,i]) $(Brownians.r[2,i]) $(Brownians.σ[i])")
    end
end

function path_check_dir(dir_name)
    path = pwd()
    dir_path = path * "/" * dir_name
    if !isdir(dir_path)
        mkdir(dir_path)
    end
    return dir_path
end

function output_physval(sys_para::systemParameters,energy::Float64,pressure::Float64,outfile::IOStream)
    println(outfile,"$(sys_para.ϕ) $(energy) $(pressure)")
end

function output_sk(sk,k_min, outfile::IOStream)
    for i = 1:length(sk)
        println(outfile,"$(i*k_min) $(sk[i])")
    end
end

function output_gr(gr,r, outfile::IOStream)
    for i = 1:length(gr)
        println(outfile,"$(r[i]) $(gr[i])")
    end
end

function output_fskt(fskt,time, outfile::IOStream)
    for i = 1:length(fskt)
        println(outfile,"$(time[i]) $(fskt[i])")
    end
end
