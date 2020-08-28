# ==============================================================================
#                           input arguments
# ==============================================================================
include("simulation_setup.jl")

function return_filename(args)
    if length(args) != 1
        println("The arguments have to be 1.
        filename should be like 'Swap_log_N1000_phi0.6_T1.0e-6_ID014.dat' ")
    end
    filename = args[1]
    return filename
end

function get_para_from_filename(filename)
    #basenameを噛ますとpathを入れた場合にもfilenameをとってくれる。
    filename = basename(filename)
    if !occursin(r"harmonic_N([\d]+)_rho([\d . ]+)_T([\d . ]+)_gamma([\d . ]+)_sigmar([\d . ]+)_([\d]+).dat",filename)
        #"harmonic_N$(simu_para.N)_rho$(simu_para.ρ)_T$(simu_para.T)_gamma$(simu_para.γ)_sigmar$(simu_para.σ_ratio)_$(str_ID).dat"

        println(filename)
        error("filename should be like 'harmonic_N1000_rho0.6_T0.1_sigmar1.4_001.dat'")
    else
        #parameter = match(r"Swap_log_N(?<number>[\d .]+)_phi(?<phi>[\d .]+)_T(?<T>[\d . - e])_ID(?<id>[\d]).dat",filename)
        parameter = match(r"harmonic_N(?<number>[\d .]+)_rho(?<rho>[\d . ]+)_T([\d . ]+)_gamma([\d . ]+)_sigmar(?<sigmar>[\d . ]+)_([\d]+).dat",filename)
        N = parse(Int64,parameter[:number])
        ρ = parse(Float64,parameter[:rho])
        σ_ratio = parse(Float64,parameter[:sigmar])
        #println(N," ",ϕ," ",T," ",id)
        return N,ρ,σ_ratio
    end
end

function open_readfile(openfile::String)::IOStream
    outfile = open(openfile,"r")
end

function read_config(filename,Brownians::Brownian,sys_para::systemParameters)
    PATH = pwd() * "/Initial_random_config/" * filename
    open(PATH) do input
        counter = 1
        while !eof(input)
            pos_x,pos_y,sigma = parse.(Float64,split(readline(input)))
            Brownians.r[1, counter] = pos_x
            Brownians.r[2, counter] = pos_y
            Brownians.σ[counter] = sigma
            if counter > sys_para.N
                println("The filesize is strange")
                break
            end
            counter +=1
        end
    end
end


function readlines_config(Brownians::Brownian,N::Int64,IOSfile::IOStream)
    for i in 1:N
        Brownians.r[1, i],Brownians.r[2, i] = parse.(Float64,split(readline(IOSfile)))
    end
    check_close(IOSfile)
end

function readlines_config_w_time(Brownians::Brownian,time::Float64,N::Int64,IOSfile::IOStream)
    for i in 1:N
        time,Brownians.r[1, i],Brownians.r[2, i] = parse.(Float64,split(readline(IOSfile)))
    end
    check_close(IOSfile)
end

function readlines_config_velocity(Brownians::Brownian,N::Int64,IOSfile_config::IOStream,IOSfile_velocity::IOStream)
    for i in 1:N
        Brownians.r[1, i],Brownians.r[2, i] = parse.(Float64,split(readline(IOSfile_config)))
        Brownians.v[1, i],Brownians.v[2, i] = parse.(Float64,split(readline(IOSfile_velocity)))
    end
    check_close(IOSfile)
end

function check_close(IOSfile)
    if eof(IOSfile)
        close(IOSfile)
    end
end

function pnum_from_filename(filename)
    file = match(r"(.*)N(?<number>\d+).*", filename)
    parse(Int64,file[:number])
end
