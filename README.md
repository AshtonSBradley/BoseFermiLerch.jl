# BoseFermi.jl

Our aim is to accurately and reliably evaluate upper incomplete Bose and Fermi integrals. We provide robust evaluation for a wide range of arguments such as occurs in thermal ultra-cold gases. 

## Now
- [x] Give definitions for upper incomplete Bose and Fermi integrals, and upper incomplete Lerch transcendent integral.
- [x] Provide reliable numerical evaluation, tested for a wide range of arguments. This is done with adaptive quadrature, and while not particularly slow, is also not optimally fast.

## Future

- [ ] implement fast evaluation using Chebychev (see e.g. [GSL](https://github.com/JuliaMath/GSL.jl)) and asymptotic expansions.  

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
While a number of asymptotic expansions are available, to provide simple and robust evaluation for a wide range of arguments, here we evaluate for $z\in \mathbb C\backslash [e^{b},\infty)$ using adaptive Gauss-Kronrod numerical quadrature provided by the [QuadGK](https://www.wikiwand.com/en/Lerch_zeta_function) package. 

### Bose and Fermi
The Bose and Fermi integrals are then evaluated via the identities:
$$g_\nu(z,\varepsilon)=z\Phi(z,\nu,1,\varepsilon)$$

$$f_\nu(z,\varepsilon)=z\Phi(-z,\nu,1,\varepsilon)$$

### Special cases 

#### $|z|<1$
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
f_s(z)=f_s(z,0)=z\Phi(-z,s,1,0)=\sum_{n=1}\frac{(-1)^nz^n}{n^s}
$$

### Riemann zeta

$$
\zeta(s)=g_s(1)=\Phi(1,s,1,0)=\sum_{n=1}^\infty\frac{1}{n^s}.
$$
