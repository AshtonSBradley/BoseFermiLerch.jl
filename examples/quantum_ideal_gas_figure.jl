using CairoMakie
using SpecialFunctions
using BoseFermiLerch

# Reproduces the three-panel quantum ideal gas comparison figure using
# BoseFermiLerch.jl as the numerical backend.

# Dimensionality parameter for a massive 3D gas.
a = 3 / 2
fcol = :dodgerblue
ccol = :mediumseagreen
bcol = :red

# FERMI GAS ========================================
zf = exp.(range(-2, 20, 201)) |> collect
zf[end] = 1e10

Tf = gamma(a + 1)^(1 / a)
Pf = Tf / (a + 1)

Tf_curve = @. fermi(zf, a)^(-1 / a)
Pf_curve = @. Tf_curve^(a + 1) * fermi(zf, a + 1)
Sf_curve = @. (a + 1) * fermi(zf, a + 1) / fermi(zf, a) - log(zf)
muf_curve = @. log(zf) * Tf_curve

Tf_curve = [Tf_curve; 0.0]
Pf_curve = [Pf_curve; Pf]
Sf_curve = [Sf_curve; 0.0]
muf_curve = [muf_curve; [Tf]]

fig = Figure(size = (960, 340), linewidth = 2)
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[1, 3])

lines!(ax1, Tf_curve, muf_curve; label = "Fermi gas", color = fcol, linewidth = 3)
text!(ax1, 0.05, 1.25; text = L"\epsilon_F")
scatter!(ax1, [Tf], [0.0]; marker = :circle, color = fcol, markersize = 12)
text!(ax1, Tf + 0.05, 0.0; text = L"T_F")

lines!(ax2, Tf_curve, Sf_curve; label = "Fermi gas", color = fcol, linewidth = 3)
lines!(ax3, Tf_curve, Pf_curve; label = "Fermi gas", color = fcol, linewidth = 3)

# CLASSICAL GAS =====================================
Tc_curve = @. zf^(-1 / a)
Pc_curve = @. Tc_curve^(a + 1) * zf
Sc_curve = @. a + 1 - log(zf)
muc_curve = @. log(zf) * Tc_curve

Tc_curve = [Tc_curve; 0.0]
Pc_curve = [Pc_curve; 0.0]
Sc_curve = [Sc_curve; -Inf]
muc_curve = [muc_curve; [0.0]]

lines!(ax1, Tc_curve, muc_curve; label = "Classical gas", color = ccol, linewidth = 3)
lines!(ax2, Tc_curve, Sc_curve; label = "Classical gas", color = ccol, linewidth = 3)
lines!(ax3, Tc_curve, Pc_curve; label = "Classical gas", color = ccol, linewidth = 3)

# BOSE GAS ==========================================
zb = range(0.01, 0.99, 99)
zb = [zb; [0.999; 0.9999; 1 - 2eps()]]

Tb_curve = @. bose(zb, a)^(-1 / a)
Pb_curve = @. Tb_curve^(a + 1) * bose(zb, a + 1)
Sb_curve = @. (a + 1) * bose(zb, a + 1) / bose(zb, a) - log(zb)
mub_curve = @. log(zb) * Tb_curve

if a > 1
    Tc = bose(1, a)^(-1 / a)
    Pc = Tc^(a + 1) * bose(1, a + 1)
    Sc = (a + 1) * bose(1, a + 1) / bose(1, a)
    T2 = range(Tc, 0, 101)
    P2 = @. T2^(a + 1) * bose(1, a + 1)
    S2 = @. T2^a * (a + 1) * bose(1, a + 1)
    mu2 = @. log(1) * T2

    Tb_curve = [Tb_curve; T2]
    Pb_curve = [Pb_curve; P2]
    Sb_curve = [Sb_curve; S2]
    mub_curve = [mub_curve; mu2]

    scatter!(ax1, [Tc], [mu2[1]]; marker = :circle, color = bcol, markersize = 12)
    text!(ax1, Tc, -0.3; text = L"T_C")
    scatter!(ax2, [Tc], [Sc]; marker = :circle, color = (bcol, 0.5), markersize = 12)
    scatter!(ax3, [Tc], [Pc]; marker = :circle, color = (bcol, 0.5), markersize = 12)
else
    Tb_curve = [Tb_curve; 0.0]
    Pb_curve = [Pb_curve; 0.0]
    Sb_curve = [Sb_curve; 0.0]
    mub_curve = [mub_curve; [0.0]]
end

lines!(ax1, Tb_curve, mub_curve; label = "Bose gas", color = bcol, linewidth = 3)
lines!(ax2, Tb_curve, Sb_curve; label = "Bose gas", color = bcol, linewidth = 3)
lines!(ax3, Tb_curve, Pb_curve; label = "Bose gas", color = bcol, linewidth = 3)

for ax in (ax1, ax2, ax3)
    xlims!(ax, 0, 1.6)
    ax.xlabel = "Temperature T/T'"
    hidedecorations!(ax; ticks = false, ticklabels = false, label = false)
end

yrange1 = 1.5 + 0.7 * a
ylims!(ax1, -0.8 * yrange1, 0.6 * yrange1)
ax1.ylabel = "Chemical potential μ/(kT')"
axislegend(ax1; position = :lb)

yrange2 = 0.2 + 1 + (1 + log(1.8)) * a
ylims!(ax2, -0.4 * yrange2, yrange2)
ax2.ylabel = "Entropy S/(Nk)"

ylims!(ax3, -0.2, 1.0)
ax3.ylabel = "Pressure P/P'"

fig
save(joinpath(@__DIR__, "quantum_ideal_gas_figure.png"), fig)
