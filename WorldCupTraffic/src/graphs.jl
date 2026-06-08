import Plots

function system_population(us, ts, ns, dxs;
        title="Population over time",
        xlabel="t",
        ylabel="Population",
        round_func= x->x)
    length(us) == length(ts) || throw(ArgumentError("time vec and pop density vecs have the same length"))

    pop = zeros(length(ts))

    for (i,u) in enumerate(us)
        offset=1
        for (n, dx) in zip(ns,dxs)
            pop[i] += sum(u[offset:offset+n-1]) * dx
            offset += n
        end
    end

    return Plots.plot(ts, round_func.(pop),
         xlabel=xlabel,
         ylabel=ylabel,
         ylims=(0.9*minimum(pop),1.1*maximum(pop)),
         title=title
    )
end

function system_energy(us, ts, ns, dxs;
        title="Energy over time",
        xlabel="t",
        ylabel="Energy",
        round_func=x->x)
    length(us) == length(ts) || throw(ArgumentError("time vec and pop density vecs have the same length"))

    energy = zeros(length(ts))

    for (i,u) in enumerate(us)
        offset=1
        for (n, dx) in zip(ns,dxs)
            energy[i] += sum(u[offset:offset+n-1].^2) * dx
            offset += n
        end
    end

    return Plots.plot(ts, round_func.(energy),
         xlabel=xlabel,
         ylabel=ylabel,
         ylims=(0.9*minimum(energy),1.1*maximum(energy)),
         title=title
    )
end

export system_population, system_energy
