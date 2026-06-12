function get_description()
    return "Test"
end

function get_model()
    μu = 0.5
    μv = 1.0
    a  = 1.5
    b  = 2.0
    pu = 0.0
    pv = 0.0

    f(u, v) = -μu * u + a * u^2 / (1.0 + v) + pu
    g(u, v) = -μv * v + b * u^2 + pv

    # Nullklina g(u,v)=0, czyli v = h(u)
    h(u) = (b * u^2 + pv) / μv

    equilibria = [
        (0.0, 0.0, "stable"),
        (0.5, 0.5, "unstable"),
        (1.0, 2.0, "DDI")
    ]

    return (
        f = f,
        g = g,
        h = h,
        equilibria = equilibria,

        u0 = -0.1,
        u1 = 2.0,
        v0 = -0.1,
        v1 = 2.5,

        Nu = 20,
        Nv = 20,
        arrow_length = 0.15,

        title = "GM basic",
        output_filename = "phase_portrait.pdf"
    )
end