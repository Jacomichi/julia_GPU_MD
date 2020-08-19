using CUDA


function ret_size(i::Int,N1::Int64,σ_ratio::Float64)
    return ifelse(i <= N1, 1.0 ,σ_ratio)
end
