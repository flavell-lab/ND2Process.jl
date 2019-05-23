function generate_MHD_str(spacing_lat, spacing_axi, size_x, size_y, size_z, filename_raw)
"ObjectType = Image
NDims = 3
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 0 1 0 0 0 1
Offset = 0 0 0
CenterOfRotation = 0 0 0
AnatomicalOrientation = RAI
ElementSpacing = $(spacing_lat) $(spacing_lat) $(spacing_axi)
DimSize = $(size_x) $(size_y) $(size_z)
ElementType = MET_USHORT
ElementDataFile = $(filename_raw)"
end

function write_MHD_str(path, MHD_str)
    open(path, "w") do f
        write(f, MHD_str)
    end
end

function write_MHD_str(path, spacing_lat, spacing_axi, size_x, size_y, size_z, filename_raw)
    MHD_str = generate_MHD_str(spacing_lat, spacing_axi, size_x, size_y, size_z, filename_raw)
    write_MHD_str(path, MHD_str)
end

function write_raw(path, array)
    open(path, "w") do f
        write(f, array)
    end
end
