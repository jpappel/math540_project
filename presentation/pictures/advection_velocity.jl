using Plots

xs = LinRange(0,1,210)
hat_vals = @. 1-max(1-abs(6 * (xs[71:end-70]-0.5)),0)
ys = [ones(70);hat_vals;ones(70)]

plt = plot(xs, ys,
     title = "Velocity",
     label = "1-max(1-|6x-0.5|,0)"
    );
savefig(plt, "advection_velocity.png");
