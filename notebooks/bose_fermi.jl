### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ f338efd3-442f-4c7a-a0da-78d897a1bbb1
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(url="https://github.com/AshtonSBradley/BoseFermi.jl.git"),
        Pkg.PackageSpec(name="Plots"),
        Pkg.PackageSpec(name="LaTeXStrings"),
		Pkg.PackageSpec(name="SpecialFunctions")
    ])
    using BoseFermi, Plots, LaTeXStrings, SpecialFunctions
end

# ╔═╡ ab942f52-3c46-11ec-181e-53e4fecd998d
md"""
# Bose-Einstien and Fermi-Dirac 
- Integrals of the BE and FD distributions
- Upper incomplete, allowing a cutoff to separation regions of phase space
- Variable precision, allowing tradeoff between speed and accuracy.
"""

# ╔═╡ 66a0b741-ff1f-4002-880c-1f992ab92ba3
md"""
Introductions
- A Primer on Quantum Fluids, Parker & Barenghi, Chapter 2
- Introduction to Statistical Physics, Kerson Huang, Chapters 9, 10, 11
"""

# ╔═╡ 988d7cf3-4477-4018-96d4-ff8035e54865
gr(size=(300,300),grid=false,legend=false)

# ╔═╡ 5d8dce92-38a1-4c60-ad9b-fbe820a4fe5d
## plot, trap, int
function new_plot(;size=(300,300))
    fs,tfs,lw = 10,8,1.5
    p = plot(legend=:left,lw=lw,size=size,frame=:box,
    foreground_color_legend = nothing,
    xtickfontsize=tfs,ytickfontsize=tfs,xguidefontsize=fs,
    yguidefontsize=fs,legendfontsize=tfs)
    return p
end

# ╔═╡ 5c485580-12c1-4045-b2ad-29b042516096
z = LinRange(0,0.999999999,600)

# ╔═╡ 679e4f75-1f7b-4e8b-b03d-bdc88fd0cf29
md"""
## Ideal Bose gas in a 3D box
For an ideal gas, the fugacity $z=e^\mu<1$. As $\mu\to0$ from below, the fugacity approaches the critical point marked by $z=1$, and the phase space density 

$$n\lambda^3=g_{3/2}(z)$$

approaches the critical value $g_{3/2}(1)=\zeta(3/2)\simeq 2.612$.
After reaching the critical density, the excited states become saturated, forcing macroscropic numbers to condense into the ground state.
"""

# ╔═╡ afc7c4c2-62ba-4c0c-b290-7ba57a2653af
begin
	p1=new_plot()
	plot!(p1,z,bose.(3/2,z))
	hline!([zeta.(3/2)],ls=:dash,legend=false)
	xlabel!(L"z")
	ylabel!(L"g_{3/2}(z)")
	ylims!(0,2.8)
	xlims!(0,1.02)
	annotate!(.3,2.49,text(String(L"\zeta(3/2)"),10))
end

# ╔═╡ a4d39e3c-ca96-4143-9cc2-f6598aad7b35
md"""
## Energy cutoff
For $T<T_c$, repulstive interactions in the Bose gas cause $\mu$ to become positive. The excited states now involve a range of energies bounded below by $\mu$. In general we can introduce a cutoff energy $\varepsilon$, such that all states with energy above $\varepsilon$ may be described by the Bose function. The requirement is now that $\mu-\varepsilon<0$"""

# ╔═╡ Cell order:
# ╟─ab942f52-3c46-11ec-181e-53e4fecd998d
# ╟─66a0b741-ff1f-4002-880c-1f992ab92ba3
# ╟─f338efd3-442f-4c7a-a0da-78d897a1bbb1
# ╠═988d7cf3-4477-4018-96d4-ff8035e54865
# ╟─5d8dce92-38a1-4c60-ad9b-fbe820a4fe5d
# ╠═5c485580-12c1-4045-b2ad-29b042516096
# ╟─679e4f75-1f7b-4e8b-b03d-bdc88fd0cf29
# ╟─afc7c4c2-62ba-4c0c-b290-7ba57a2653af
# ╟─a4d39e3c-ca96-4143-9cc2-f6598aad7b35
