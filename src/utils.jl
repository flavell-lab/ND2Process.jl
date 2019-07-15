function create_dir(dirpath::String)
    if !isdir(dirpath)
        mkdir(dirpath)
    end
end

function maxprj(array; dims)
    dropdims(maximum(array, dims=dims), dims=dims)
end

function encode_movie(input, output; fps=10)
    input = `-i $input`
    options = `-hide_banner -loglevel panic -y -framerate $fps`
    codec = `-c:v libx264 -crf 16 -pix_fmt yuv420p -preset veryslow`
    run(`ffmpeg $options $input $codec $output`)
    nothing
end
