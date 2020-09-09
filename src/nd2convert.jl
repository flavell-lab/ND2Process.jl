"""
    nd2_to_mhd(path_nd2, path_save,
        spacing_lat, spacing_axi, generate_MIP::Bool;
        θ=nothing, x_crop::Union{Nothing, UnitRange{Int64}}=nothing,
        y_crop::Union{Nothing, UnitRange{Int64}}=nothing,
        z_crop::Union{Nothing, UnitRange{Int64}}=nothing, chs::Array{Int}=[1],
        MHD_dir_name="MHD", MIP_dir_name="MIP")

Saves nd2 into MHD files after rotating and cropping. Rotation is skipped if
θ is set to `nothing`.

Arguments
---------
* `path_nd2`: path of .nd2 file to use
* `path_save`: path of .h5 file to save
* `spacing_lat`: lateral spacing (for logging)
* `spacing_axi`: axial spacing (for logging)
* `generate_MIP`: if true, save MIP in as preview
* `θ`: yaw angle (lateral rotation, radian). nothing if no rotation
* `x_crop`: x range to use. Full range if nothing
* `y_crop`: y range to use. Full range if nothing
* `z_crop`: z range to use. Full range if nothing
* `chs`: ch to use
* `MHD_dir_name`: name of the subfolder to save MHD files
* `MIP_dir_name`: name of the subfolder to save MIP files
* `n_bin`: number of rounds to bin. e.g. `n_bin=2` results in 4x4 binning
* `z_range`: number of frames per z-stack, if using a continuous timestream data series
"""
function nd2_to_mhd(path_nd2, path_save,
    spacing_lat, spacing_axi, generate_MIP::Bool;
    θ=nothing, x_crop::Union{Nothing, UnitRange{Int64}}=nothing,
    y_crop::Union{Nothing, UnitRange{Int64}}=nothing,
    z_crop::Union{Nothing, UnitRange{Int64}}=nothing, chs::Array{Int}=[1],
    MHD_dir_name="MHD", MIP_dir_name="MIP", n_bin=nothing, z_range=nothing)

    mhd_paths = []
    x_size, y_size, z_size, t_size, c_size = nd2dim(path_nd2)

    if !isnothing(z_range)
        z_size = z_range
        t_size = t_size ÷ z_range
    end
    if !isnothing(n_bin)
        x_size = floor(Int, x_size / (2 ^ n_bin))
        y_size = floor(Int, y_size / (2 ^ n_bin))
    end

    # directories
    f_basename = splitext(basename(path_nd2))[1]

    path_dir_MHD = joinpath(path_save, MHD_dir_name)
    path_dir_MIP = joinpath(path_save, MIP_dir_name)

    create_dir(path_save)
    create_dir(path_dir_MHD)
    generate_MIP && create_dir(path_dir_MIP)

    x_size_save = Int(0)
    y_size_save = Int(0)
    z_size_save = Int(0)

    if z_crop == nothing
        z_size_save = z_size
        z_crop = 1:z_size
    else
        z_size_save = length(z_crop)
    end
    if x_crop == nothing
        x_size_save = x_size
        x_crop = 1:x_size
    else
        x_size_save = length(x_crop)
    end
    if y_crop == nothing
        y_size_save = y_size
        y_crop = 1:y_size
    else
        y_size_save = length(y_crop)
    end


    img_ = zeros(UInt16, x_size, y_size)
    vol_ = zeros(UInt16, x_size_save, y_size_save, z_size_save)

    @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
        @showprogress for t_ = 1:t_size
        for c_ = chs
            for (n_, z_) = enumerate(z_crop)
                # load
                if isnothing(z_range)
                    img_ = Float64.(transpose(images.get_frame_2D(c=c_-1,
                        t=t_-1, z=z_-1)))
                else
                    img_ = Float64.(transpose(images.get_frame_2D(c=c_-1,
                        t=0, z=z_range*(t_-1)+z_-1)))
                end
                # binning
                if !isnothing(n_bin)
                    img_ = bin_img(img_, n_bin)
                end

                # rotate, crop, convert to UInt16
                if isnothing(θ)
                    vol_[:,:,n_] = round.(UInt16, img_[x_crop, y_crop])
                else
                    vol_[:,:,n_] = round.(UInt16,
                        rotate_img(img_, θ)[x_crop, y_crop])
                end
            end

            save_basename = f_basename * "_t$(lpad(t_, 4, "0"))_ch$(c_)"
            path_file_MHD = joinpath(path_dir_MHD, save_basename * ".mhd")
            path_file_raw = joinpath(path_dir_MHD, save_basename * ".raw")

            # save MHD
            write_raw(path_file_raw, vol_)
            write_MHD_spec(path_file_MHD, spacing_lat, spacing_axi,
                    x_size_save, y_size_save, z_size_save,
                        save_basename * ".raw")

            # save MIP
            if generate_MIP
                path_file_MIP = joinpath(path_dir_MIP, save_basename * ".png")
                imsave(path_file_MIP, dropdims(maximum(vol_, dims=3), dims=3),
                    cmap="gray")
            end
        end

        end
    end
end

"""
    nd2_to_h5(path_nd2, path_save, spacing_lat, spacing_axi; θ=nothing,
        x_crop::Union{Nothing, UnitRange{Int64}}=nothing,
        y_crop::Union{Nothing, UnitRange{Int64}}=nothing,
        z_crop::Union{Nothing, UnitRange{Int64}}=nothing, chs::Array{Int}=[1])

Saves nd2 into HDF5 file after rotating and cropping. Rotation is skipped if
θ is set to `nothing`. Note: indexing is 1 based. Array axis is in the following
order: [x, y, z, t, c]. HDF5 file is chunked with:
(x size, y size, z size, 1, 1, 1)

Arguments
---------
* `path_nd2`: path of .nd2 file to use
* `path_save`: path of .h5 file to save
* `spacing_lat`: lateral spacing (for logging)
* `spacing_axi`: axial spacing (for logging)
* `θ`: yaw angle (lateral rotation, rauab). nothing if no rotation
* `x_crop`: x range to use. Full range if nothing
* `y_crop`: y range to use. Full range if nothing
* `z_crop`: z range to use. Full range if nothing
* `chs`: ch to use
* `n_bin`: number of rounds to bin. e.g. `n_bin=2` results in 4x4 binning
* `z_range`: number of frames per z-stack, if using a continuous timestream data series
"""
function nd2_to_h5(path_nd2, path_save, spacing_lat, spacing_axi; θ=nothing,
    x_crop::Union{Nothing, UnitRange{Int64}}=nothing,
    y_crop::Union{Nothing, UnitRange{Int64}}=nothing,
    z_crop::Union{Nothing, UnitRange{Int64}}=nothing, chs::Array{Int}=[1],
    n_bin=nothing, z_range=nothing)

    if splitext(path_save)[2] != ".h5"
        error("path_save must end with .h5")
    end

    x_size, y_size, z_size, t_size, c_size = nd2dim(path_nd2)
    
    if !isnothing(z_range)
        z_size = z_range
        t_size = t_size ÷ z_range
    end
    if !isnothing(n_bin)
        x_size = floor(Int, x_size / (2 ^ n_bin))
        y_size = floor(Int, y_size / (2 ^ n_bin))
    end

    x_size_save = Int(0)
    y_size_save = Int(0)
    z_size_save = Int(0)

    if z_crop == nothing
        z_size_save = z_size
        z_crop = 1:z_size
    else
        z_size_save = length(z_crop)
    end
    if x_crop == nothing
        x_size_save = x_size
        x_crop = 1:x_size
    else
        x_size_save = length(x_crop)
    end
    if y_crop == nothing
        y_size_save = y_size
        y_crop = 1:y_size
    else
        y_size_save = length(y_crop)
    end


    img_ = zeros(UInt16, x_size, y_size)

    @pywith py_nd2reader.ND2Reader(path_nd2) as images begin

        h5open(path_save, "w") do f
            dset = d_create(f, "data", datatype(UInt16),
            dataspace(x_size_save, y_size_save, z_size_save, t_size,
                length(chs)), "chunk", (x_size_save, y_size_save, 1, 1, 1))

            @showprogress for t_ = 1:t_size
                for (n_c, c_) = enumerate(chs)
                    for (n_z, z_) = enumerate(z_crop)
                        # load
                        if isnothing(z_range)
                            img_ = Float64.(transpose(images.get_frame_2D(c=c_-1,
                                t=t_-1, z=z_-1)))
                        else
                            img_ = Float64.(transpose(images.get_frame_2D(c=c_-1,
                                t=0, z=z_range*(t_-1)+z_-1)))
                        end

                        # binning
                        if !isnothing(n_bin)
                            img_ = bin_img(img_, n_bin)
                        end

                        # rotate, crop, convert to UInt16, and save
                        if isnothing(θ)
                            dset[:, :, n_z, t_, n_c] = round.(UInt16,
                            img_[x_crop, y_crop])
                        else
                            dset[:, :, n_z, t_, n_c] = round.(UInt16,
                                rotate_img(img_, θ)[x_crop, y_crop])
                        end
                    end # for
                end # for
            end # for
        end # h5open

        if isnothing(θ) θ = 0.0 end
        dict_attr = Dict("x_crop"=>[x_crop.start, x_crop.stop],
            "y_crop"=>[y_crop.start, y_crop.stop],
            "z_crop"=>[z_crop.start, z_crop.stop], "θ"=>θ,
            "spacing_lat"=>spacing_lat, "spacing_axi"=>spacing_axi)
        h5writeattr(path_save, "data", dict_attr)
    end # pywith
end

"""
    function write_nd2_preview(path_nd2; prjdim=3, chs=[1], z_crop=:drop_first)

Saves maximum intensity projection (MIP) of nd2 file and make movies of the
time series.

Arguments
---------
* `path_nd2`: path of .nd2 file to use
* `prjdim`: MIP projection dimension. Default: 3 (across z slices)
* `z_crop`: z range to use. Full range if nothing. `:drop_first` drops 1st frame
* `chs`: ch to use
* `dir_save`: directory to save MIP images and movies
* `n_bin`: number of rounds to bin. e.g. `n_bin=2` results in 4x4 binning
* `z_range`: number of frames per z-stack, if using a continuous timestream data series
"""
function write_nd2_preview(path_nd2; prjdim=3, chs=[1], z_crop=:drop_first,
    dir_save=nothing, n_bin=nothing, z_range=nothing)
    x_size, y_size, z_size, t_size, c_size = nd2dim(path_nd2)

    if !isnothing(n_bin)
        x_size = floor(Int, x_size / (2 ^ n_bin))
        y_size = floor(Int, y_size / (2 ^ n_bin))
    end

    if isnothing(dir_save)
        dir_save = dirname(path_nd2)
    end

    dir_MIP = joinpath(dir_save, "MIP_original")
    dir_movie = joinpath(dir_save, "movie_original")
    f_basename = splitext(basename(path_nd2))[1]

    create_dir(dir_MIP)
    create_dir(dir_movie)

    z_size_save = Int(0)
    if !isnothing(z_range)
        z_size = z_range
        t_size = t_size ÷ z_range
    end
    if z_crop == nothing
        z_size_save = z_size
        z_crop = 1:z_size
    elseif z_crop == :drop_first
        z_size_save = z_size - 1
        z_crop = 2:z_size
    else
        z_size_save = length(z_crop)
    end

    vol_ = zeros(Float32, x_size, y_size, z_size_save)
    img_MIP_ = zeros(Float32, x_size, y_size)
    n_leading_zero = 4

    # png generation
    @pywith py_nd2reader.ND2Reader(path_nd2) as images begin
        @showprogress for t_ = 1:t_size
            for (n_c, c_) = enumerate(chs)
                for (n_z, z_) = enumerate(z_crop)
                    # load
		    if isnothing(z_range)
                        img_ = Float32.(transpose(images.get_frame_2D(c=c_-1,
                            t=t_-1, z=z_-1)))
                    else
                        img_ = Float32.(transpose(images.get_frame_2D(c=c_-1,
                            t=0, z=z_range*(t_-1)+z_-1)))
                    end
                    if !isnothing(n_bin)
                        img_ = bin_img(img_, n_bin)
                    end

                    vol_[:,:,n_z] = img_
                end #z

                img_MIP_ = maxprj(vol_, dims=prjdim)

                name_png_ = f_basename * "_c" * lpad(string(c_), 2, "0") *
                    "_t" * lpad(string(t_), n_leading_zero, "0") * ".png"
                path_png_ = joinpath(dir_MIP, name_png_)
                imsave(path_png_, img_MIP_, cmap="gray", vmax=1000)
            end # c
        end # t
    end # pywith

    # movie generation
    try
        for c_ = chs
            for fps_ = [30,60,120,240]
                encode_movie(joinpath(dir_MIP, f_basename * "_c" *
                    lpad(string(c_), 2, "0") * "_t%0$(n_leading_zero)d.png"),
                    joinpath(dir_movie, f_basename *
                        "_original_$(fps_)fps.mp4"), fps=fps_);
            end
        end
    catch
        error("Cannot generate movie for $path_nd2")
    end
end
