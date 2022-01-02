# calculate reservoir parameters
# μ, ecut, for given N, T, V(r)
# using Pkg
# pkg"activate ."

using BoseFermi
using BenchmarkTools, Test

≈(a,b) = isapprox(a,b,rtol=1e-4)

## Lerchphi test 
z,s,a,b =0.1,1.5,1.,0.

## bose
@btime bose(z,s,b,rtol=1e-9)
@btime z*lerch(z,s,a,b,rtol=1e-9)

# fermi
@btime fermi(z,s,b,rtol=1e-9)
@btime z*lerch(-z,s,a,b,rtol=1e-9)

z,s,a,b =0.5,1.0,1.,0.
z*lerch(z,s,a,b)
bose(z,s,b)

z,s,a,b = 1.0,2.0,1.,0.
@btime z*lerch(z,s,a,b)
@btime bose(z,s,b)

## simple timing
using Plots, LaTeXStrings
using SpecialFunctions
using HCubature
using Roots

    @btime bose(1,3/2,0.)
    @btime zeta(3/2)
    @btime bose(.5,3,.5)
    exp(.5)
    zeta(3)
    bose(.99,0.1,.1)

    # note only getting 5 digits for log(2)
    @test bose(1,3)==bose(1,3,0)
    @test bose(1,5)==zeta(5)
    @test bose(.5,1)≈log(2)

    μ = 12. # in oscillator units
    β = 1/(6μ) #temp in units of μ
    ecut = 3*μ

    @btime bose(exp(β*μ),3/2,β*ecut)

    exp(β*μ),exp(.1)
    V0(x,y,z) = .5*x^2 + .5*y^2 + .5*z^2
    ΔV(x,y) = 10*exp(-.5(x^2 + y^2))

## Total potential
    V(x,y,z) = V0(x,y,z) + ΔV(x,y)

    # Effective potential (Hartree-Fock)
    Veff(x,y,z,μ) = V(x,y,z) < μ ? 2μ - V(x,y,z) : V(x,y,z)

    # examine the effective potential
    xcut = sqrt(2*ecut)
    xi = 0.
    xf = .7xcut
    Nx = 500
    g = 0.1
    μ = 11
    x = LinRange(xi,xf,Nx) |> Vector
    Vmin = minimum(V.(x,0.,0.))
    ntf_plot(x,μ) = V(x,0.,0.) < μ ? (μ - V(x,0.,0.))/g + Vmin/g : Vmin/g
    δμ = μ - Vmin # Thomas-Fermi particle density
    ntf(x,μ) = V(x...) < μ ? (μ - V(x...))/g : 0.0

## The TF denstiy must be defined with
# effective μ relative to the trap minimum, i.e. μ -> μ_rel = μ - V_min

    plot(x,Veff.(x,0.,0.,δμ),label=L"V_\textrm{eff}(r)",grid=false)
    plot!(x,V.(x,0.,0.),label=L"V(r)")
    plot!(x,V0.(x,0.,0.),line=:dash,label=L"V_0(r)")
    plot!(x,ΔV.(x,0.),line=:dash,label=L"\delta V(r)")
    plot!(x,g*ntf_plot.(x,δμ),line=:dashdot,label=L"gn_{\tiny\textrm{TF}}(r)")
    plot!(x,x./x*μ,label=L"\mu",legend=:bottomright)
    plot!(x,x./x*δμ,label=L"\delta\mu",legend=:bottomright)

    # define cutoffs
    Kc(x,y,z) = sqrt(2*maximum([ecut - V(x,y,z),0.0]))
    Kceff(x,y,z,μ) = sqrt(2*maximum([ecut - Veff(x,y,z,μ),0.0]))


    xcut = sqrt(2*ecut)
    xi = 0.
    xf = 2.3*xcut
    Nx = 500
    x = LinRange(xi,xf,Nx) |> Vector

    #slice along x for I-region density
    plot(x,bose.(exp.(β*(μ.-V.(x,0.,0.))),3/2),label="no cutoff")
    plot!(x,bose.(exp.(β*(μ.-V.(x,0.,0.))),3/2,β*ecut),label=L"\epsilon_{cut}\textrm{ (wrong cutoff)}")
    plot!(x,bose.(exp.(β*(μ.-V.(x,0.,0.))),3/2,β*Kc.(x,0.,0.).^2/2),label=L"\hbar^2K_{cut}(\mathbf{r})^2/2m",ylabel=L"n(\mathbf{r})\lambda_{dB}^3",xlabel=L"r",title="Reservoir particle density")
    #plot(x,bose.(3/2,exp.(β*(μ-Veff.(x,0.,0.,μ)))),label="Hartree-Fock")
    # annotate(L"g_{3/2}(e^{\beta(\mu-V(\mathbf{r}))})",[1.4,2.3])
    # annotate(L"g_{3/2}(e^{\beta(\mu-V(\mathbf{r}))},\frac{\beta \hbar^2}{2m}K_{cut}(\mathbf{r})^2)",[.5,1.4])
    # annotate(L"g_{3/2}(e^{\beta(\mu-V(\mathbf{r}))},\beta \epsilon_{cut})",[.7,.7])
    # xlim(0,19);ylim(0,2.5)


    # Calculate chemical potential and cutoff energy
    # We use the Hartree-Fock approach, but leave the potential entirely general: density of states not required!

## finite domain test
    xmin = (-20.,-20.,-20.)
    xmax = (20.,20.,20.)
    nth(x,μ) = bose(exp(β*(μ-Veff(x...,μ))),3/2,β*Kc.(x...).^2/2)
    numeric,err = hcubature(x->nth(x,μ), xmin, xmax,rtol=1e-3)

    Nth(μ) = hcubature(x->nth(x,μ), xmin, xmax,rtol=1e-3)[1]

    # @code_warntype Nth(10.)

    @time Nth(10);
    Nth(16)
    g = 0.1

## TODO
xmin = (-6.,-6.,-6.)
xmax = (6.,6.,6.)
ntf(x,μ) = V(x...) < μ ? (μ - V(x...))/g : 0.0
Ntf(μ) = hcubature(x->ntf(x,μ), xmin, xmax,rtol=1e-3)[1]

Ntf(10)
Ntf(16)

Ntot(μ) = Ntf(μ) + Nth(μ)
Ntot(10)
Ntot(16)
NT = 3e4

## TODO Tidy up Thomas-Fermi limits, and bose limits for numerical integration

μp = LinRange(10,16,10) |> Vector
p = plot(grid=false)

plot!(μp,Ntot.(μp),marker="o",label=L"N_{tot}")
plot!(μp,μp./μp*NT,label=L"N_{target}")
plot!(μp,Ntf.(μp),marker="o",label=L"N_{TF}")
plot!(μp,Nth.(μp),marker="o",label=L"N_{th}")
plot!(p,xlabel=L"\mu",ylabel=L"N(\mu)",legend=:bottomright)

μ0 = find_zero(x->Ntot(x)-NT, (12,16))

plot(μp,Nth.(μp),marker="o",label=L"N_{th}")
