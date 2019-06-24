function nd2read(path_nd2; ch=1, t=1)
    x_size, y_size, c_size, t_size, z_size = nd2dim(path_nd2)

    stack_ = (length(t) == 1) ? zeros(UInt16, z_size, x_size, y_size) : zeros(
        UInt16, z_size, x_size, y_size, length(t))

    if length(t) > 1 # multiple time point
        @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
            for (n_, t_) = enumerate(t)
                for z_ = 1:z_size
                    stack_[z_,:,:,n_] = transpose(images.get_frame_2D(c=ch-1,
                        t=t_-1, z=z_-1))
                end
            end
        end
    else # 1 time point
        @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
            for z_ = 1:z_size
                stack_[z_,:,:] = transpose(images.get_frame_2D(c=ch-1, t=t-1,
                    z=z_-1))
            end
        end
    end

    stack_
end

function nd2dim(path_nd2)
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
    x_size, y_size, c_size, t_size, z_size = nd2dim(path_nd2)

    t_list = [1, round(Int, t_size / 2), t_size]
    stack_ = zeros(UInt16, z_size, x_size, y_size, 3)
    @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
        for (n_, t_) = enumerate(t_list)
            for z_ = 1:z_size
                stack_[z_,:,:,n_] = transpose(images.get_frame_2D(c=ch-1, t=t_-1, z=z_-1))
            end
        end
    end

    stack_MIP = Float64.(dropdims(maximum(stack_, dims=1), dims=1))

    for i = 1:3
        subplot(1,3,i)
        imshow(stack_MIP[:,:,i])
        title("t=$(t_list[i]+1)")
    end
    tight_layout()

    return_data && return stack_
end

function nd2preview(stack::Array; θ, x_crop=nothing, y_crop=nothing, z_crop=nothing)
    @assert size(stack, 4) == 3

    if z_crop == nothing
        z_crop = 1:size(stack, 1)
    end
    if y_crop == nothing
        y_crop = 1:size(stack, 2)
    end
    if x_crop == nothing
        x_crop = 1:size(stack, 3)
    end


    stack_MIP = Float64.(dropdims(maximum(stack[z_crop,:,:,:], dims=1), dims=1))
    stack_MIP_proc = zeros(eltype(stack_MIP), length(x_crop), length(y_crop), 3)

    for i = 1:3
        stack_MIP_proc[:,:,i] = rotate_img(stack_MIP[:,:,i], θ)[x_crop, y_crop]
    end

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
