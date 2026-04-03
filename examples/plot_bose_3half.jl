using CairoMakie
using BoseFermiLerch

z = range(0.0, 1.0 - 1e-6; length = 400)
g = bose.(z, 3 / 2)

fig = Figure(size = (700, 420))
ax = Axis(
    fig[1, 1];
    xlabel = "z",
    ylabel = L"g_{3/2}(z)",
    title = L"Bose\ function\ g_{3/2}(z) = \mathrm{Li}_{3/2}(z)",
)

lines!(ax, z, g; color = :firebrick, linewidth = 3)
xlims!(ax, 0, 1)

fig
save(joinpath(@__DIR__, "bose_3half.png"), fig)
