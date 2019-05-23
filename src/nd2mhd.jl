function nd2_dim(path_nd2)
    @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
        @assert eltype(images.get_frame_2D(c=0,t=0,z=0)) == UInt16
        x_size, y_size, c_size, t_size, z_size = [images.sizes[k] for k = ["x", "y", "c", "t", "z"]]
        return (x_size, y_size, c_size, t_size, z_size)
    end
end

function rotate_img(img, θ)
    tfm = recenter(RotMatrix(θ), center(img))
    ImageTransformations.warp(img, tfm)
end

function nd2preview(path_nd2; ch=1, return_data=false)
    x_size, y_size, c_size, t_size, z_size = nd2_dim(path_nd2)

    t_list = [0, round(Int, t_size / 2), t_size-1]
    stack_ = zeros(UInt16, z_size, x_size, y_size, 3)
    @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
        for (n_, t_) = enumerate(t_list)
            for z_ = 1:z_size
                stack_[z_,:,:,n_] = transpose(images.get_frame_2D(c=ch-1, t=t_, z=z_-1))
            end
        end
    end

    stack_MIP = Float64.(dropdims(maximum(stack_, dims=1), dims=1))

    figure()
    for i = 1:3
        subplot(1,3,i)
        imshow(stack_MIP[:,:,i])
        title("t=$(t_list[i]+1)")
    end
    tight_layout()

    return_data && return stack_
end

function nd2preview(stack::Array; θ, x_crop, y_crop, z_crop=nothing)
    @assert size(stack, 4) == 3

    if z_crop == nothing
        z_crop = 1:size(stack, 1)
    end

    stack_MIP = Float64.(dropdims(maximum(stack[z_crop,:,:,:], dims=1), dims=1))
    stack_MIP_proc = zeros(eltype(stack_MIP), length(x_crop), length(y_crop), 3)

    for i = 1:3
        stack_MIP_proc[:,:,i] = rotate_img(stack_MIP[:,:,i], θ)[x_crop, y_crop]
    end

    figure()
    for i = 1:3
        subplot(1,3,i)
        imshow(stack_MIP_proc[:,:,i])
    end
    tight_layout()

    if true in isnan.(stack_MIP_proc)
        @warn "transformed data contains NaN due to rotation."
    end

    println(size(stack_MIP_proc, 1) ./ [2,4,8,16])
    println(size(stack_MIP_proc, 2) ./ [2,4,8,16])
end

function nd2_to_mhd(path_nd2, path_save,
    spacing_lat, spacing_axi, generate_MIP::Bool
    θ, x_crop, y_crop, z_crop=nothing, chs::Array{Int}=[1],
    MHD_dir_name="MHD", MIP_dir_name="MIP")

    mhd_paths = []
    x_size, y_size, c_size, t_size, z_size = nd2_dim(path_nd2)

    # directories
    f_basename = splitext(basename(path_nd2))[1]

    path_dir_MHD = joinpath(path_save, MHD_dir_name)
    path_dir_MIP = joinpath(path_save, MIP_dir_name)

    create_dir(path_save)
    create_dir(path_dir_MHD)
    generate_MIP && create_dir(path_dir_MIP)

    x_size_save = length(x_crop)
    y_size_save = length(y_crop)
    z_size_save = Int(0)

    if z_crop == nothing
        z_size_save = z_size
        z_crop = 1:z_size
    else
        z_size_save = length(z_crop)
    end

    img_ = zeros(UInt16, x_size, y_size)
    vol_ = zeros(UInt16, x_size_save, y_size_save, z_size_save)

    @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
        @showprogress for t_ = 1:t_size
        for c_ = chs
            for (n_, z_) = enumerate(z_crop)
                # load
                img_ = Float64.(transpose(images.get_frame_2D(c=c_-1, t=t_-1, z=z_-1)))
                # rotate and crop
                vol_[:,:,n_] = round.(UInt16, rotate_img(img_, θ)[x_crop, y_crop])
            end

            save_basename = f_basename * "_t$(lpad(t_, 4, "0"))_ch$(c_)"
            path_file_MHD = joinpath(path_dir_MHD, save_basename * ".mhd")
            path_file_raw = joinpath(path_dir_MHD, save_basename * ".raw")

            # save MHD
            write_raw(path_file_raw, permutedims(vol_, [2,1,3]))
            write_MHD_str(path_file_MHD, spacing_lat, spacing_axi,
                    y_size_save, x_size_save, z_size_save, save_basename * ".raw")

            # save MIP
            if generate_MIP
                path_file_MIP = joinpath(path_dir_MIP, save_basename * ".png")
                imsave(path_file_MIP, dropdims(maximum(vol_, dims=3), dims=3), cmap="gray")
            end
        end

        end
    end
end
