using Plots

xs = LinRange(0, 2*pi, 200)
ys = sin.(xs)

plt = plot(xs, ys,
     title = "Test Plot",
     label = "sin(x)"
    );
savefig(plt, "test.png");
