using Sockets, Serialization

function get_snapshot(host::AbstractString = "127.0.0.1", port::Int = 20000)
    sock = connect(host, port)
    snap = deserialize(sock)
    close(sock)
    return snap
end

get_snapshot()
# => powinno zwrócić (:hello, [0.123, ...], 1.23)