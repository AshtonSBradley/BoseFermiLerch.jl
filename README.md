# BoseFermiLerch.jl

Upper incomplete Bose and Fermi integrals. Robust evaluation for a wide range of arguments such as occurs in thermal ultra-cold gases, based on the [Lerch transcendent](https://en.wikipedia.org/wiki/Lerch_zeta_function) and its upper incomplete integral extension. 

![](/examples/quantum_ideal_gas_figure.png)

## Evaluation strategy

This package aims for a robust default evaluation strategy across both real and
complex arguments. For the complete Lerch transcendent (`b = 0`), it uses a
hybrid map: the convergent gamma-series in the easy `|z| < 1` regime where it
is broadly valid and typically cheapest, a direct real-axis integral on the
negative real line to stabilize large-fugacity Fermi-Dirac evaluations, and a
contour-integral backend in the spirit of
[Computing the Lerch transcendent](https://fredrikj.net/blog/2022/02/computing-the-lerch-transcendent/).
The contour method acts as the general fallback, including for the complete
positive-real near-`z = 1` regime where the series becomes expensive. For the
upper incomplete case (`b > 0`), the implementation first reuses the same
series path when `|z| < 1` is still favourable, and otherwise combines the
complete evaluation with a tail correction computed by adaptive quadrature.

In practice, this gives a good tradeoff between speed and reliable evaluation,
especially near difficult regions such as branch cuts, large fugacity, and
complex arguments. It should also serve as a strong fallback baseline for future
optimisation strategies.

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

While a number of asymptotic expansions are available, to provide simple and
robust evaluation for a wide range of arguments, this package evaluates for
$z\in \mathbb C\backslash [e^{b},\infty)$. The practical evaluation map is:

- for broad `|z| < 1` regimes, use the convergent series below
- for complete negative-real arguments, use direct adaptive quadrature on the real axis
- otherwise, use the contour backend, and for `b > 0` subtract the lower-tail correction

This keeps the easy cases fast while preserving a robust fallback for difficult
near-cut and complex inputs.

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

## Speed regression baseline

The repository also includes a small performance regression harness under `perf/`.
A recent strict run produced the following median timings:

| Case | Arguments `(z, s, a, b)` | Median time (ms) | Budget (ms) |
| --- | --- | ---: | ---: |
| `complete_well_damped` | `(0.5, 2, 1, 0)` | `0.002` | `0.400` |
| `complete_near_one` | `(0.999999, 2, 1, 0)` | `0.009` | `1.500` |
| `incomplete_small_b` | `(0.5, 2, 1, 0.05)` | `0.001` | `1.500` |
| `incomplete_near_one_small_b` | `(0.999, 2, 1, 0.05)` | `0.002` | `1.500` |
| `incomplete_just_above_one_small_b` | `(1.01, 2, 1, 0.05)` | `0.002` | `1.500` |
| `incomplete_large_b` | `(0.5, 2, 1, 1)` | `0.001` | `1.500` |
| `real_cut_window` | `(1.5, 2, 1, 1)` | `0.002` | `1.500` |
| `fermi_complete_equiv` | `(-0.8, 2, 1, 0)` | `0.002` | `0.600` |
| `complex_complete` | `(0.6 + 0.2im, 2, 1, 0)` | `0.001` | `2.500` |
| `complex_incomplete` | `(0.6 + 0.2im, 2, 1, 0.5)` | `0.002` | `3.000` |

To rerun the strict check:

```bash
BFL_ENFORCE_SPEED=1 julia --project=perf -e 'include("perf/speed_regression.jl")'
```
