using GLMakie

function TwoNodeArcAnim(
    sol,
    ns,
    Ls;
    max_frames=300,
    scale=0.5,
    filename="graph_transport.mp4"
)

    n_frames = length(sol.t)

    if n_frames > max_frames
        step = ceil(Int, n_frames / max_frames)
        indices = 1:step:n_frames
    else
        indices = 1:n_frames
    end

    times = sol.t[indices]

    θuv = collect(LinRange(0, π, ns[1]))
    θvu = collect(LinRange(π, 2π, ns[2]))

    Puv0 = @view sol.u[indices[1]][1:ns[1]]
    Pvu0 = @view sol.u[indices[1]][ns[1]+1:ns[1]+ns[2]]

    uv_points = Observable(
        Point2f[
            (
                (1 + scale*Puv0[i]) * cos(θuv[i]),
                (1 + scale*Puv0[i]) * sin(θuv[i])
            )
            for i in eachindex(θuv)
        ]
    )

    vu_points = Observable(
        Point2f[
            (
                (1 + scale*Pvu0[i]) * cos(θvu[i]),
                (1 + scale*Pvu0[i]) * sin(θvu[i])
            )
            for i in eachindex(θvu)
        ]
    )

    fig = Figure(size=(800,800))

    ax = Axis(
        fig[1,1],
        aspect=DataAspect(),
        title="t=$(round(times[1], digits=2))"
    )

    hidedecorations!(ax)

    # Reference graph
    lines!(
        ax,
        cos.(LinRange(0, π, 200)),
        sin.(LinRange(0, π, 200)),
        linestyle=:dash
    )

    lines!(
        ax,
        cos.(LinRange(π, 2π, 200)),
        sin.(LinRange(π, 2π, 200)),
        linestyle=:dash
    )

    # Nodes
    scatter!(ax, [-1, 1], [0, 0], markersize=20)

    lines!(ax, uv_points, linewidth=4, label="u → v")
    lines!(ax, vu_points, linewidth=4, label="v → u")

    limits!(ax, -2, 2, -2, 2)

    record(fig, filename, eachindex(indices)) do k

        Puv = @view sol.u[indices[k]][1:ns[1]]
        Pvu = @view sol.u[indices[k]][ns[1]+1:ns[1]+ns[2]]

        uv_points[] = Point2f[
            (
                (1 + scale*Puv[i]) * cos(θuv[i]),
                (1 + scale*Puv[i]) * sin(θuv[i])
            )
            for i in eachindex(θuv)
        ]

        vu_points[] = Point2f[
            (
                (1 + scale*Pvu[i]) * cos(θvu[i]),
                (1 + scale*Pvu[i]) * sin(θvu[i])
            )
            for i in eachindex(θvu)
        ]

        ax.title = "t = $(round(times[k], digits=2))"
    end

    return fig
end

export TwoNodeArcAnim
