include("md_struct.jl")

function return_ID(args)
    if length(args) != 1
        println("The arguments have to be 1.
        ID::Int64")
    end
    ID = parse(Int64,args[1])
    return ID
end
