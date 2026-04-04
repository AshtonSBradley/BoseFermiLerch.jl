using CairoMakie
using BoseFermiLerch
using LaTeXStrings
using SpecialFunctions

function styled_axis(figpos)
    return Axis(
        figpos;
        xlabel = L"z",
        ylabel = L"n\lambda^3",
        xlabelsize = 20,
        ylabelsize = 20,
        xticklabelsize = 16,
        yticklabelsize = 16,
        topspinevisible = true,
        rightspinevisible = true,
    )
end

z = range(0.0, 0.999999999; length = 600)
ec = 0.1
zc = exp(ec)
zcut = range(0.0, 0.999999999; length = 600)

fig = Figure(size = (820, 360))

ax1 = styled_axis(fig[1, 1])
lines!(ax1, z, bose.(z, 3 / 2); color = :dodgerblue, linewidth = 3)
hlines!(ax1, [zeta(3 / 2)]; color = :black, linestyle = :dash, linewidth = 1.5)
xlims!(ax1, 0, 1.02)
ylims!(ax1, 0, 2.8)
text!(ax1, 0.30, 2.34; text = L"\zeta(3/2)", fontsize = 18)
text!(ax1, 0.55, 1.00; text = L"g_{3/2}(z)", fontsize = 18)
scatter!(ax1, [0.999999], [zeta(3 / 2)]; color = :red, markersize = 10)

ax2 = styled_axis(fig[1, 2])
lines!(ax2, zcut, bose.(zcut, 3 / 2, ec); color = :dodgerblue, linewidth = 3)
vlines!(ax2, [zc]; color = :black, linestyle = :dash, linewidth = 1.5)
xlims!(ax2, 0, 1.02)
ylims!(ax2, 0, 2.8)
text!(ax2, 0.55, 1.10; text = L"g_{3/2}(z,\varepsilon)", fontsize = 18)
text!(ax2, 0.78, 3.20; text = L"e^{\varepsilon}", fontsize = 18)

fig
save(joinpath(@__DIR__, "bose_3half.png"), fig)
