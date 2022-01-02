# BoseFermi.jl
Evaluate upper incomplete Bose and Fermi integrals. The aim is to provide robust accurate evaluation for a wide range of arguments such as occurs in thermal ultra-cold gases. 

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

### Incomplete Lerch transcendent
The Bose and Fermi integrals are defined in terms of the upper incomplete [Lerch transcendent](https://www.wikiwand.com/en/Lerch_zeta_function) defined here as
$$
\Phi(z,s,a,b)\equiv\frac{1}{\Gamma(s)}\int_b^\infty \frac{t^{s-1}e^{-at}}{1-ze^{-t}}dt
$$
While a number of asymptotic expansions are available, to provide simple and robust evaluation, here we simply evaluate for $z\in \mathbb C\backslash [e^{b},\infty)$ using adaptive Gauss-Kronrod numerical quadrature provided by the [QuadGK](https://www.wikiwand.com/en/Lerch_zeta_function) package. 

Then
$$g_\nu(z,\varepsilon)=z\Phi(z,\nu,1,\varepsilon)$$

$$f_\nu(z,\varepsilon)=z\Phi(-z,\nu,1,\varepsilon)$$

## Examples
### Bose

### Fermi
