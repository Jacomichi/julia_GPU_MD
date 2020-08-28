# ==============================================================================
#                         Equation of motion
# ==============================================================================

include("./simulation_setup.jl")
include("./simulation_observe.jl")
include("./simulation_output.jl")
include("./simulation_input.jl")

function main()
    filename = return_filename(ARGS)
    N,ρ,σ_ratio = get_para_from_filename(filename)

    sys_para = systemParameters(N,ρ,σ_ratio)
    Brownians = Brownian(sys_para.N,sys_para.σ_ratio)

    #setting analyze
    k_min = 2.0*π/sys_para.L
    k_max = 2.0*π

    dr = 0.1
    r = setting_gr(sys_para,dr)

    IOSdata = open_readfile(filename)

    counter = 0
    while !eof(IOSdata)
        readlines_config(Brownians,sys_para.N,IOSdata)
        if counter == 0
            sk = calc_sk(Brownians,sys_para,k_min,k_max)
            gr = calc_one_gr(r,dr,Brownians,sys_para)
        else
            sk += calc_sk(Brownians,sys_para,k_min,k_max)
            gr += calc_one_gr(r,dr,Brownians,sys_para)
        end
        counter += 1
    end

    sk = sk./counter
    sk_filename = "sk_" * basename(filename)
    sk_outfile = open_writefile(sk_filename)
    output_sk(sk,k_min, sk_outfile)

    gr = gr./counter
    gr_filename = "gr_" * basename(filename)
    gr_outfile = open_writefile(gr_filename)
    output_gr(gr,r,gr_outfile)
end

main()
