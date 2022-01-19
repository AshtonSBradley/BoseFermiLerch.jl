# BoseFermi.jl

Our aim is to accurately and reliably evaluate upper incomplete Bose and Fermi integrals. We provide robust evaluation for a wide range of arguments such as occurs in thermal ultra-cold gases. 

## Now
- [x] Give definitions for upper incomplete Bose and Fermi integrals, and upper incomplete Lerch transcendent integral.
- [x] Provide reliable numerical evaluation for a wide range of arguments. This is done with adaptive quadrature, and while not particularly slow, is also not optimally fast 

## Future

- [ ] implement fast evaluation using Chebychev and asymptotic expansions.  

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

### Incomplete Lerch transcendent
The Bose and Fermi integrals may be defined in terms of the more general upper incomplete [Lerch transcendent](https://www.wikiwand.com/en/Lerch_zeta_function) defined here as
$$
\Phi(z,s,a,b)\equiv\frac{1}{\Gamma(s)}\int_b^\infty \frac{t^{s-1}e^{-at}}{1-ze^{-t}}dt
$$
While a number of asymptotic expansions are available, to provide simple and robust evaluation for a wide range of arguments, here we evaluate for $z\in \mathbb C\backslash [e^{b},\infty)$ using adaptive Gauss-Kronrod numerical quadrature provided by the [QuadGK](https://www.wikiwand.com/en/Lerch_zeta_function) package. 

The Bose and Fermi integrals are then evaluated via the identities:
$$g_\nu(z,\varepsilon)=z\Phi(z,\nu,1,\varepsilon)$$

$$f_\nu(z,\varepsilon)=z\Phi(-z,\nu,1,\varepsilon)$$

## Examples
### Bose

### Fermi

### Lerch 