using GLMakie
import LinearAlgebra: norm

function directed_arc(p1::Point2f, p2::Point2f, sign, npts)

    d = p2 - p1
    L = norm(d)

    midpoint = (p1 + p2) / 2

    normal = Point2f(-d[2], d[1])
    normal /= norm(normal)

    center = midpoint + sign * 0.35 * L * normal

    θ1 = atan(p1[2] - center[2], p1[1] - center[1])
    θ2 = atan(p2[2] - center[2], p2[1] - center[1])

    if sign > 0
        θ2 < θ1 && (θ2 += 2π)
    else
        θ2 > θ1 && (θ2 -= 2π)
    end

    θ = range(θ1, θ2, length=npts)

    r = norm(p1 - center)

    Point2f[
        center + r * Point2f(cos(t), sin(t))
        for t in θ
    ]
end

function curve_normals(curve)

    N = length(curve)

    normals = Vector{Point2f}(undef, N)

    for i in 1:N

        tangent =
            if i == 1
                curve[2] - curve[1]
            elseif i == N
                curve[N] - curve[N-1]
            else
                curve[i+1] - curve[i-1]
            end

        n = Point2f(-tangent[2], tangent[1])
        normals[i] = n / norm(n)
    end

    normals
end

function FourNodeArcAnim(
    sol,
    Nx;
    scale=0.25,
    max_frames=300,
    filename="graph_transport.mp4"
)

    #
    # Node locations
    #

    node_pos = Dict(
        1 => Point2f(0.0, 0.0),
        2 => Point2f(2.0, 0.0),
        3 => Point2f(1.0, 1.75),
        4 => Point2f(0.0, -2.0)
    )

    #
    # Directed edges in same order as state vector
    #

    edges = [
        (1,2),
        (2,1),

        (1,3),
        (3,1),

        (1,4),
        (4,1),

        (2,3),
        (3,2)
    ]

    #
    # Build reference arcs
    #

    reference_paths = [
        # straight line
        [Point2f(t*node_pos[r][1]+(1-t)*node_pos[s][1], t*node_pos[r][2]+(1-t)*node_pos[s][2]) for t in LinRange(0,1,Nx)]
        for (k, (s,r)) in enumerate(edges) 
    ]
    
    #=
    reference_paths = [
        directed_arc(
            node_pos[s],
            node_pos[t],
            isodd(k) ? 1 : -1,
            Nx
        )
        for (k,(s,t)) in enumerate(edges)
    ]
    =#

    reference_normals =
        curve_normals.(reference_paths)

    #
    # Frame selection
    #

    n_frames = length(sol.t)

    if n_frames > max_frames
        step = ceil(Int, n_frames/max_frames)
        indices = 1:step:n_frames
    else
        indices = 1:n_frames
    end

    times = sol.t[indices]

    #
    # Initial geometry
    #

    initial_state = sol.u[first(indices)]

    edge_obs = Observable[]

    for edge_idx in 1:length(edges)

        density =
            initial_state[
                (edge_idx-1)*Nx+1 : edge_idx*Nx
            ]

        pts = Point2f[
            reference_paths[edge_idx][i] +
            scale *
            density[i] *
            reference_normals[edge_idx][i]
            for i in 1:Nx
        ]

        push!(edge_obs, Observable(pts))
    end

    #
    # Figure
    #

    fig = Figure(size=(900,900))

    ax = Axis(
        fig[1,1],
        aspect=DataAspect(),
        title="t = $(round(times[1], digits=2))"
    )

    hidedecorations!(ax)
    hidespines!(ax)

    #
    # Reference graph
    #

    for path in reference_paths

        lines!(
            ax,
            first.(path),
            last.(path),
            linestyle=:dash,
            linewidth=2
        )
    end

    #
    # Animated transport curves
    #

    for obs in edge_obs

        lines!(
            ax,
            obs,
            linewidth=4
        )
    end

    #
    # Nodes
    #

    scatter!(
        ax,
        [node_pos[i][1] for i in 1:4],
        [node_pos[i][2] for i in 1:4],
        markersize=20
    )

    text!(ax, 0.0,  0.15, text="1")
    text!(ax, 2.0,  0.15, text="2")
    text!(ax, 1.0,  1.95, text="3")
    text!(ax, 0.0, -2.2, text="4")

    limits!(ax, -1.5, 3.0, -3.0, 3.0)

    #
    # Animation
    #

    record(fig, filename, eachindex(indices)) do frame

        state = sol.u[indices[frame]]

        for edge_idx in 1:length(edges)

            density =
                state[
                    (edge_idx-1)*Nx+1 : edge_idx*Nx
                ]

            edge_obs[edge_idx][] = Point2f[
                reference_paths[edge_idx][i] +
                scale *
                density[i] *
                reference_normals[edge_idx][i]
                for i in 1:Nx
            ]
        end

        ax.title =
            "t = $(round(times[frame], digits=2))"
    end

    return fig
end

export FourNodeArcAnim
