const py_nd2reader = PyNULL()

function __init__()
    copy!(py_nd2reader, pyimport_conda("nd2reader", "nd2reader",
        "conda-forge"))
end
