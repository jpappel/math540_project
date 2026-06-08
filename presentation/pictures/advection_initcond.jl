using Plots

xs = LinRange(0,1,210)
ys = max.(cos.(4*pi*xs),0)

plt = plot(xs, ys,
     title = "Initial Conditions",
     label = "max(cos(4πx),0)"
    );
savefig(plt, "advection_initcond.png");
