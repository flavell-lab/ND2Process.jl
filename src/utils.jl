function create_dir(dirpath::String)
    if !isdir(dirpath)
        mkdir(dirpath)
    end
end

function maxprj(array; dims)
    dropdims(maximum(array, dims=dims), dims=dims)
end
