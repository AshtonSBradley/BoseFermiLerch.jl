import Pkg
Pkg.activate(temp=true)
Pkg.add(["BoseFermiLerch", "CairoMakie", "QuadGK"])

using CairoMakie
using QuadGK
using SpecialFunctions
using BoseFermiLerch

function big_bose_ref(z, s; rtol = big"1e-30")
    setprecision(256) do
        zb = BigFloat(z)
        sb = BigFloat(s)
        gb = gamma(sb)
        integrand(t) = exp(-t) * t^(sb - 1) * zb / ((1 - zb * exp(-t)) * gb)
        quadgk(integrand, big"0.0", big"1.0", big"4.0", big"16.0", big"64.0", Inf; rtol = rtol)[1]
    end
end

zs = 1 .- 10.0 .^ range(-5, -2; length = 60)
orders = [1.5, 2.5]
colors = [:steelblue, :darkorange]

fig = Figure(size = (1100, 460))
ax_err = Axis(fig[1, 1], xlabel = "1 - z", ylabel = "relative error", xscale = log10, yscale = log10, title = "Bose near z = 1: asymptotic patch accuracy")
ax_time = Axis(fig[1, 2], xlabel = "1 - z", ylabel = "time (ns)", xscale = log10, yscale = log10, title = "Bose near z = 1: asymptotic vs contour")

for (s, color) in zip(orders, colors)
    asymp_err = Float64[]
    public_ns = Float64[]
    contour_ns = Float64[]

    for z in zs
        ref = Float64(big_bose_ref(z, s))
        value = bose(z, s; rtol = 1e-12)
        push!(asymp_err, abs(value - ref) / abs(ref))
        push!(public_ns, @elapsed(bose(z, s; rtol = 1e-12)) * 1e9)
        push!(contour_ns, @elapsed(z * BoseFermiLerch.lerch_complete_contour(z, s, 1.0; rtol = 1e-12)) * 1e9)
    end

    lines!(ax_err, 1 .- zs, asymp_err; color = color, label = "s = $(s)")
    lines!(ax_time, 1 .- zs, public_ns; color = color, label = "public s = $(s)")
    lines!(ax_time, 1 .- zs, contour_ns; color = color, linestyle = :dash)
end

axislegend(ax_err, position = :lb)
axislegend(ax_time, position = :rt)
save(joinpath(@__DIR__, "bose_near_one.png"), fig)
