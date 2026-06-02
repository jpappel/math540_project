using DifferentialEquations
using Plots

function coupled_transport!(dP, P, params, t)
    ns, dxs, velocities, boundary_conds, attractions = params

    # u -> uv edge
    # v -> vu edge
    Pu = @view P[1:ns[1]]
    Pv = @view P[ns[1]+1:ns[1]+ns[2]]
    # time derrivatives
    dPu = @view dP[1:ns[1]]
    dPv = @view dP[ns[1]+1:ns[1]+ns[2]]

    cu, cv = velocities
    dxu, dxv = dxs
    # division is expensive!
    inv_dxu = 1/dxu
    inv_dxv = 1/dxv

    total_P = sum(P)
    Au = attractions[1](t)
    Av = attractions[2](t)

    @inbounds for i in 2:ns[1]
        x = (i-1)*dxu
        # PERF: can half the number of computations by swapping cur with prev at end of loop
        #       this trades vectorization for space function evaluation of cu
        F_cur = cu(x,t,total_P,Au,Av) * Pu[i]
        F_prev = cu(x-dxu,t,total_P,Au,Av) * Pu[i-1]
        dPu[i] = -(F_cur - F_prev) * inv_dxu
    end

    @inbounds for i in 2:ns[2]
        x = (i-1)*dxv
        # PERF: see previous loop
        F_cur = cv(x,t,total_P,Av,Au) * Pv[i]
        F_prev  = cv(x-dxv,t,total_P,Av,Au) * Pv[i-1]
        dPv[i] = -(F_cur - F_prev) * inv_dxv
    end

    if isnothing(boundary_conds)
        dPu[1] = -(cu(0,t,total_P,Au,Av) * Pu[1] - cv(ns[2]*dxv,t,total_P,Av,Au) * Pv[end]) * inv_dxu
        dPv[1] = -(cv(0,t,total_P,Av,Au) * Pv[1] - cu(ns[1]*dxu,t,total_P,Au,Av) * Pu[end]) * inv_dxv
    else
        dPu[1] = -(cu(0,t,total_P,Au,Av) * Pu[1] - cv(ns[2]*dxv,t,total_P,Av,Au) * Pv[end]) * inv_dxu + boundary_conds[1](t)
        dPv[1] = -(cv(0,t,total_P,Av,Au) * Pv[1] - cu(ns[1]*dxu,t,total_P,Au,Av) * Pu[end]) * inv_dxv + boundary_conds[1](t)
    end
end

function constant_v(x,t,P, Astart, Aend)
    return 1
end

function TwoNodeSystem(;
        n1::Int=10,
        n2::Int=10,
        L1::T=1.0,
        L2::T=1.0,
        init_P1::Vector{Float64}=ones(n1)/n1,
        init_P2::Vector{Float64}=ones(n2)/n2,
        velocity_1::Function=constant_v,
        velocity_2::Function=constant_v,
        boundary_conds::Union{Nothing, Tuple}=nothing,
        attraction_1::Function=(t)->1.0,
        attraction_2::Function=(t)->1.0,
        tspan::Tuple{T,T}=(0.0, 1.0)
    ) where {T<:Real}
    n1 == length(init_P1) || throw(ArgumentError("n1 must match length of init_P1"))
    n2 == length(init_P2) || throw(ArgumentError("n2 must match length of init_P2"))
    L1 > 0 || throw(ArgumentError("length L1 must be positive"))
    L2 > 0 || throw(ArgumentError("length L2 must be positive"))

    ns = [n1, n2]
    Ls = [L1, L2]
    dxs = [L1/n1, L2/n2]
    init_conds = [init_P1, init_P2]
    velocities = (velocity_1, velocity_2)
    attractions = (attraction_1, attraction_2)

    params = (ns, dxs, velocities, boundary_conds, attractions)
    P0 = vcat(init_P1, init_P2)

    return ODEProblem(coupled_transport!, P0, tspan, params)
end

function TwoNodePlot(ns, Ls, Pu, Pv;title="Two Node System")
    thetau = collect(LinRange(0, pi, ns[1]))
    thetav = collect(LinRange(pi, 2*pi, ns[2]))

    plt = plot(
        thetau, Pu,
        label = "uv",
        proj = :polar,
        title = title,
    )
    plot!(plt,
        thetav, Pv,
        label = "vu"
    )

    return plt
end

function TwoNodeAnim(sol, ns, Ls; max_frames=300)
    n_frames = length(sol.t)

    if n_frames > max_frames
        step = ceil(Int, n_frames / max_frames)
        indicies = 1:step:n_frames
    else
        indicies = 1:n_frames
    end
    times = sol.t[indicies]

    all_Puv = [@view sol.u[k][1:ns[1]] for k in indicies]
    all_Pvu = [@view sol.u[k][ns[1]+1:ns[1]+ns[2]] for k in indicies]

    thetau = collect(LinRange(0, pi, ns[1]))
    thetav = collect(LinRange(pi, 2*pi, ns[2]))

    maxP = maximum(vcat(
        [maximum(p) for p in all_Puv],
        [maximum(p) for p in all_Pvu]
    ))

    # base plot to avoid expensive reallocs
    plt = plot(
        proj=:polar,
        legend=:topright,
        ylims=(0, maxP * 1.05),
        title="t=$(round(sol.t[1], digits=2))"
    )
    plot!(plt, thetau, all_Puv[1], label="uv")
    plot!(plt, thetav, all_Pvu[1], label="vu")

    anim = @animate for k in 1:length(times)
        # Update data (faster than recreating)
        plt.series_list[1][:y] = all_Puv[k]
        plt.series_list[2][:y] = all_Pvu[k]
        title!(plt, "t=$(round(times[k], digits=2))")
        plt
    end

    return anim
end

export TwoNodeSystem, TwoNodeAnim, TwoNodePlot
