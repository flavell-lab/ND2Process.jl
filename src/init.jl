const py_nd2reader = PyNULL()

function __init__()
    copy!(py_nd2reader, pyimport("nd2reader"))
end
