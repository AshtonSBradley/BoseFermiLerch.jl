module BoseFermi

import SpecialFunctions: zeta, gamma, gamma_inc

export bose, fermi, lerch 

gamma_inc_Q(p,x) = gamma_inc(p,x)[2]

function lerch_term(z,s,a,b,k)
    return z^k/(k+a)^s*gamma_inc_Q(s,(k+a)*b)
end

function lerch(z,s,a=1.0,b=0.0;rtol=1e-8,atol=rtol)
    @assert 0.0 <= s
    @assert 0.0 <= b
    @assert 0.0 < a
    (z == 1 && b == 0) && return zeta(s)
    (z == 0.5 && s == 1) && return 2*log(2)
    k = 0
    Sk = lerch_term(z,s,a,b,k)
    Skp1 = Sk + lerch_term(z,s,a,b,k+1)
    while !isapprox(Sk,Skp1,rtol=rtol,atol=atol)
        k += 1
        Sk = Skp1
        Skp1 += lerch_term(z,s,a,b,k+1)
    end
    return Sk
end

"""
`bose(z,ν,y)`

Evaluates the incomplete Bose-Einstein function

```math
g_\\nu(z,y)=\\frac{1}{\\Gamma(\\nu)}\\int_{y}^\\infty dx\\frac{x^{\\nu-1}}{z^{-1}e^{x}-1} = \\frac{1}{\\Gamma(\\nu)}\\sum_{k=1}^\\infty \\frac{z^k}{k^{\\nu}} \\Gamma(\\nu,ky).
```

When \$y=0\$, reduces to the regular Bose-Einstein function, equivalent to the polylogarithm:

```math
g_\\nu(z,0)= g_\\nu(z)\\equiv \\text{Li}_\\nu(z).
```

For \$z=1\$, \$y=0\$, reduces to the Reimann-zeta function
```math
g_\\nu(1,0)=\\zeta(\\nu)=\\sum_{k=1}^\\infty \\frac{1}{k^\\nu}.
```
This implementation requires the normalized incomplete gamma function.
"""
function bose(z,s,b=0;rtol=1e-9,atol=rtol) 
    return z*lerch(z,s,1.0,b,rtol=rtol,atol=atol)
end

"""
`fermi(z,ν,y)`

Evaluates the incomplete Fermi-Dirac function

```math
f_\\nu(z,y)=\\frac{1}{\\Gamma(\\nu)}\\int_{y}^\\infty dx\\frac{x^{\\nu-1}}{z^{-1}e^{x}+1} = \\frac{1}{\\Gamma(\\nu)}\\sum_{k=1}^\\infty \\frac{(-z)^k}{k^{\\nu}} \\Gamma(\\nu,ky).
```

When \$y=0\$, reduces to the regular Bose-Einstein function, equivalent to the polylogarithm:

```math
f_\\nu(z,0)= f_\\nu(z)\\equiv \\text{Li}_\\nu(-z).
```

For \$z=1\$, \$y=0\$, reduces to the Reimann-zeta function
```math
f_\\nu(1,0)=\\zeta(\\nu)=\\sum_{k=1}^\\infty \\frac{1}{k^\\nu}.
```
This implementation requires the normalized incomplete gamma function.
"""
function fermi(z,s,b=0;rtol=1e-9,atol=rtol)
    return z*lerch(-z,s,1.0,b,rtol=rtol,atol=atol)
end 

end