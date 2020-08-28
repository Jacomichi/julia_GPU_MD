include("md_struct.jl")

function create_filename_harmonic_config(simu_para::Simu_Param;ID::Int64)::String
    str_ID = lpad(ID,3,'0')
    filename = "harmonic_N$(simu_para.N)_rho$(simu_para.ρ)_T$(simu_para.T)_gamma$(simu_para.γ)_sigmar$(simu_para.σ_ratio)_$(str_ID).dat"
    return filename
end

function create_filename_config(words::String,simu_para::Simu_Param;ID::Int64)::String
    str_ID = lpad(ID,3,'0')
    filename = "$(words)_N$(simu_para.N)_rho$(simu_para.ρ)_T$(simu_para.T)_gamma$(simu_para.γ)_sigmar$(simu_para.σ_ratio)_$(str_ID).dat"
    return filename
end


function path_check_dir(dir_name)
    path = pwd()
    dir_path = path * "/" * dir_name
    if !isdir(dir_path)
        mkdir(dir_path)
    end
    return dir_path
end

function create_file_harmonic_config(simu_para::Simu_Param,id::Int64=0)::String
    path_dir = path_check_dir("harmonic_config")
    filename = create_filename_harmonic_config(simu_para,ID=id)
    config_openfile = path_dir * "/" * filename
    return config_openfile
end

function create_file_config(dir_name::String,filename_header::String,simu_para::Simu_Param,id::Int64=0)::String
    path_dir = path_check_dir(dir_name)
    filename = create_filename_config(filename_header,simu_para,ID=id)
    config_openfile = path_dir * "/" * filename
    return config_openfile
end

function open_writefile(openfile::String)::IOStream
    outfile = open(openfile,"w")
end

function output_position(pos, outfile::IOStream)
    for i = 1:size(pos)[2]
        println(outfile,"$(pos[1,i]) $(pos[2,i])")
    end
end
