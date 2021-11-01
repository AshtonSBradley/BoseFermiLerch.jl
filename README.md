# BoseFermi.jl
Evaluate upper incomplete Bose and Fermi integrals. The strategy here is to provide variable precision for distributed field applications, such as occurs in thermal ultra-cold gases.

## Definitions
### Bose-Einstein functions
$$
g_\nu(z,\varepsilon)=\frac{1}{\Gamma(\nu)}\int_\varepsilon^\infty dx\frac{x^{\nu-1}}{e^xz^{-1}-1}
$$

### Fermi-Dirac functions
$$
f_\nu(z,\varepsilon)=\frac{1}{\Gamma(\nu)}\int_\varepsilon^\infty dx\frac{x^{\nu-1}}{e^xz^{-1}+1}
$$

where the fugacity $z\equiv e^{\mu}$, and the lower integration is restricted to $\varepsilon>\mu$.

## Examples
### Bose

### Fermi
