# BoseFermiLerch.jl

Upper incomplete Bose and Fermi integrals. Robust evaluation for a wide range of arguments such as occurs in thermal ultra-cold gases, based on the [Lerch transcendent](https://en.wikipedia.org/wiki/Lerch_zeta_function) and its upper incomplete integral extension.

## Now
- [x] Give definitions for upper incomplete Bose and Fermi integrals, and upper incomplete Lerch transcendent.
- [x] Reliable numerical evaluation using a convergent series for `|z| < 1` and adaptive quadrature elsewhere.

## Future

- [ ] fast evaluation using Chebychev (see e.g. [GSL](https://github.com/JuliaMath/GSL.jl)) and asymptotic expansions.  

## Definitions
### Bose-Einstein integrals
$$
g_\nu(z,\varepsilon)=\frac{1}{\Gamma(\nu)}\int_\varepsilon^\infty dx\frac{x^{\nu-1}}{e^xz^{-1}-1}
$$

### Fermi-Dirac integrals
$$
f_\nu(z,\varepsilon)=\frac{1}{\Gamma(\nu)}\int_\varepsilon^\infty dx\frac{x^{\nu-1}}{e^xz^{-1}+1}
$$

where the fugacity $z\equiv e^{\mu}$, and the lower integration is restricted to $\varepsilon>\mu$.

## Incomplete Lerch integral
The Bose and Fermi integrals may be defined in terms of the more general upper incomplete [Lerch integral](https://www.wikiwand.com/en/Lerch_zeta_function) defined here as

$$
\Phi(z,s,a,b)\equiv\frac{1}{\Gamma(s)}\int_b^\infty \frac{t^{s-1}e^{-at}}{1-ze^{-t}}dt
$$

While a number of asymptotic expansions are available, to provide simple and robust evaluation for a wide range of arguments, this package evaluates for $z\in \mathbb C\backslash [e^{b},\infty)$. For $|z|<1$ it uses the convergent series below; otherwise it uses adaptive Gauss-Kronrod numerical quadrature from [QuadGK](https://juliamath.github.io/QuadGK.jl/stable/).

### Bose and Fermi
The Bose and Fermi integrals are then evaluated via the identities:
$$g_\nu(z,\varepsilon)=z\Phi(z,\nu,1,\varepsilon)$$

$$f_\nu(z,\varepsilon)=z\Phi(-z,\nu,1,\varepsilon)$$

### Special cases 

#### |z|<1
In this case, Lerch integral can be written as the convergent series

$$ 
\Phi(z,s,a,b)\equiv\frac{1}{\Gamma(s)}\int_b^\infty \frac{t^{s-1}e^{-at}}{1-ze^{-t}}dt=\frac{1}{\Gamma(s)}\sum_{n=0}^\infty \frac{z^n}{(a+n)^s}\Gamma(s,b(a+n))
$$

### Lerch transcendent
For $b=0$ this reduces to the [Lerch transcendent](https://en.wikipedia.org/wiki/Lerch_zeta_function)

$$
\Phi(z,s,a,0)=\Phi(z,s,a)\equiv\sum_{n=0}^\infty \frac{z^n}{(a+n)^s}.
$$

### Incomplete Bose function
$$
g_s(z,b)=z\Phi(z,s,1,b)=\frac{1}{\Gamma(s)}\sum_{n=1}\frac{z^n}{n^s}\Gamma(s,bn)
$$

### Bose function
$$
g_s(z)=g_s(z,0)=z\Phi(z,s,1,0)=\sum_{n=1}\frac{z^n}{n^s}
$$

### Fermi function 
$$
f_s(z)=f_s(z,0)=z\Phi(-z,s,1,0)=\sum_{n=1}^\infty\frac{(-1)^{n-1}z^n}{n^s}=-\operatorname{Li}_s(-z)
$$

### Special values

$$
\zeta(s)=g_s(1)=\Phi(1,s,1,0)=\sum_{n=1}^\infty\frac{1}{n^s}.
$$

$$
f_s(1)=\eta(s)=(1-2^{1-s})\zeta(s).
$$

## Domain notes

- `bose(z, s, b)` and `fermi(z, s, b)` require `s > 0` and `b >= 0`.
- `lerch(z, s, a, b)` evaluates the principal branch away from the real branch cut `z in [exp(b), Inf)`.
- `fermi(z, s, b)` inherits that restriction through `-z`, so it excludes real `z <= -exp(b)`.

## Example

An example plot of `bose(z, 3/2)` for `z in [0, 1)` is included in
`examples/plot_bose_3half.jl`.

One way to run it is:

```julia
import Pkg
Pkg.activate(temp=true)
Pkg.add(["BoseFermiLerch", "CairoMakie"])
include("examples/plot_bose_3half.jl")
```

The script writes `examples/bose_3half.png`.

A larger example reproducing the three-panel quantum ideal gas comparison figure is
included in `examples/quantum_ideal_gas_figure.jl`. It writes
`examples/quantum_ideal_gas_figure.png`.
