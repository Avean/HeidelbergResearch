using Sockets, Serialization

last_snapshot = Ref((:hello, rand(3), 1.23))  # cokolwiek seryjnego

function snapshot_server!(port::Int = 20000)
    server = listen(port)
    @async begin
        while true
            sock = accept(server)
            snap = deepcopy(last_snapshot[])
            serialize(sock, snap)
            flush(sock)
            close(sock)
        end
    end
    return server
end

snapshot_server!()