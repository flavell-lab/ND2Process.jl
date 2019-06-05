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
