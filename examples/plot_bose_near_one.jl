import Pkg
Pkg.activate(temp = true)
Pkg.add(["BoseFermiLerch", "CairoMakie", "QuadGK", "SpecialFunctions", "Chairmarks"])

using BoseFermiLerch
using CairoMakie
using Chairmarks
using QuadGK
using SpecialFunctions
using Statistics

const FIGDIR = joinpath(@__DIR__, "figures")
mkpath(FIGDIR)

const NONINTEGER_ORDERS = (0.5, 1.5, 2.5)
const INTEGER_ORDERS = (1.0, 2.0, 3.0)
const NEAR_ONE_Z = (0.99, 0.999, 0.9999, 0.99999)
const SWITCH_Z = (0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
const MONOTONE_Z = (0.9, 0.99, 0.999, 0.9999)

function big_bose_ref(z, s; rtol = big"1e-30", prec = 256)
    @assert 0 < z < 1
    @assert s > 0
    setprecision(prec) do
        zb = BigFloat(z)
        sb = BigFloat(s)
        gb = gamma(sb)
        integrand(t) = zb * exp(-t) * t^(sb - 1) / ((1 - zb * exp(-t)) * gb)
        quadgk(
            integrand,
            big"0.0", big"1.0", big"4.0", big"16.0", big"64.0", Inf;
            rtol = rtol,
        )[1]
    end
end

relerr(x, y) = abs(x - y) / abs(y)

function public_bose(z, s)
    return bose(z, s; rtol = 1e-12)
end

function bose_reference(z, s)
    if s == 1.0
        return -log1p(-z)
    end
    return Float64(big_bose_ref(z, s))
end

function chair_median_time_ns(f; seconds = 0.15)
    bench = @be ($f)() seconds = seconds
    return median(getproperty.(bench.samples, :time)) * 1e9
end

function savefig_pdf_png(fig, stem::AbstractString)
    save(joinpath(FIGDIR, stem * ".pdf"), fig)
    save(joinpath(FIGDIR, stem * ".png"), fig)
end

function figure_test_results()
    fig = Figure(size = (1200, 860))

    ax_nonint = Axis(
        fig[1, 1],
        xlabel = L"1-z",
        ylabel = "relative error",
        xscale = log10,
        yscale = log10,
        title = "Noninteger near-one accuracy",
    )
    hlines!(ax_nonint, [1e-10], color = :black, linestyle = :dash, label = "test tolerance")

    for s in NONINTEGER_ORDERS
        errs = Float64[]
        for z in NEAR_ONE_Z
            ypkg = public_bose(z, s)
            yref = bose_reference(z, s)
            push!(errs, relerr(ypkg, yref))
        end
        lines!(ax_nonint, 1 .- collect(NEAR_ONE_Z), errs, label = "s = $(s)")
        scatter!(ax_nonint, 1 .- collect(NEAR_ONE_Z), errs)
    end
    axislegend(ax_nonint, position = :rb)

    ax_int = Axis(
        fig[1, 2],
        xlabel = L"1-z",
        ylabel = "relative error",
        xscale = log10,
        yscale = log10,
        title = "Integer near-one accuracy",
    )
    hlines!(ax_int, [1e-10], color = :black, linestyle = :dash, label = "s = 2,3 tolerance")
    hlines!(ax_int, [1e-13], color = :gray40, linestyle = :dot, label = "s = 1 tolerance")

    for s in INTEGER_ORDERS
        errs = Float64[]
        for z in NEAR_ONE_Z
            ypkg = public_bose(z, s)
            yref = bose_reference(z, s)
            push!(errs, relerr(ypkg, yref))
        end
        lines!(ax_int, 1 .- collect(NEAR_ONE_Z), errs, label = "s = $(Int(s))")
        scatter!(ax_int, 1 .- collect(NEAR_ONE_Z), errs)
    end
    axislegend(ax_int, position = :rb)

    ax_switch = Axis(
        fig[2, 1],
        xlabel = "z",
        ylabel = "relative error",
        yscale = log10,
        title = "Switching-region continuity",
    )
    hlines!(ax_switch, [5e-10], color = :black, linestyle = :dash, label = "test tolerance")
    if isdefined(BoseFermiLerch, :Z1_REAL_ASYMP_MU_CUTOFF)
        zswitch = exp(-getfield(BoseFermiLerch, :Z1_REAL_ASYMP_MU_CUTOFF))
        vlines!(ax_switch, [zswitch], color = :gray50, linestyle = :dot, label = "activation threshold")
    end

    for s in (0.5, 1.5)
        errs = Float64[]
        for z in SWITCH_Z
            ypkg = public_bose(z, s)
            yref = bose_reference(z, s)
            push!(errs, relerr(ypkg, yref))
        end
        lines!(ax_switch, collect(SWITCH_Z), errs, label = "s = $(s)")
        scatter!(ax_switch, collect(SWITCH_Z), errs)
    end
    axislegend(ax_switch, position = :rb)

    ax_mono = Axis(
        fig[2, 2],
        xlabel = "z",
        ylabel = L"g_s(z)",
        title = "Monotonicity sanity check",
    )
    for s in (NONINTEGER_ORDERS..., INTEGER_ORDERS...)
        vals = [public_bose(z, s) for z in MONOTONE_Z]
        lines!(ax_mono, collect(MONOTONE_Z), vals, label = "s = $(s)")
        scatter!(ax_mono, collect(MONOTONE_Z), vals)
    end
    axislegend(ax_mono, position = :lt)

    Label(
        fig[0, :],
        "Near-z=1 Bose public-test results",
        fontsize = 24,
        font = :bold,
        padding = (0, 0, 12, 0),
    )

    savefig_pdf_png(fig, "fig_bose_near_one_test_results")
end

function figure_runtime_public()
    fig = Figure(size = (900, 600))
    ax = Axis(
        fig[1, 1],
        xlabel = L"1-z",
        ylabel = "median runtime (ns/eval)",
        xscale = log10,
        yscale = log10,
        title = "Runtime of public bose(z,s) on the near-one test grid",
    )

    for s in NONINTEGER_ORDERS
        tmed = Float64[]
        for z in NEAR_ONE_Z
            push!(tmed, chair_median_time_ns(() -> public_bose(z, s)))
        end
        lines!(ax, 1 .- collect(NEAR_ONE_Z), tmed, label = "s = $(s)")
        scatter!(ax, 1 .- collect(NEAR_ONE_Z), tmed)
    end

    axislegend(ax, position = :lb)
    savefig_pdf_png(fig, "fig_bose_near_one_runtime")
end

figure_test_results()
figure_runtime_public()
