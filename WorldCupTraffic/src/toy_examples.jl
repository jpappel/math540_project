using DifferentialEquations

function two_node(cs,N, init_conds, bondary_conds, tspan=(0.0,1.0))
    dx = 1/(N-1)
    xs = 0:dx:1

    function coupled_transport!()
    end

    prob = ODEProblem(init_conds, tspan, p)
end


# Initial condition
u0 = sin.(2π .* x)

# Boundary condition
g(t) = 0.0

function advection!(du, u, p, t)
    c, dx = p

    # Left boundary (Dirichlet)
    u_left = g(t)

    # Upwind scheme for c > 0
    du[1] = -c * (u[1] - u_left)/dx

    for i in 2:length(u)
        du[i] = -c * (u[i] - u[i-1])/dx
    end
end

tspan = (0.0, 1.0)

prob = ODEProblem(advection!, u0, tspan, (c, dx))

sol = solve(prob, Tsit5(), saveat=0.01)
