module BoseFermi

import SpecialFunctions: zeta, gamma 
using QuadGK 

export bose, fermi, lerch 

lerch_int(t,z,s,a) = 1/gamma(s)t^(s-1)*exp(-a*t)/(1-z*exp(-t))

lerch(z,s,a,b;rtol=1e-8) = quadgk(t->lerch_int(t,z,s,a),b,Inf,rtol=rtol)[1]


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
function bose(z,s,b=0;rtol=1e-9) 
    @assert 0 <= s
    @assert 0 <= b
    if z == one(z) && b == 0
        return zeta(s)
    elseif z == 0.5 && s == 1
        return log(2)
    else 
        return z*lerch(z,s,1.0,b,rtol=rtol)
    end
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
function fermi(z,s,b=0;rtol=1e-9)
    @assert 0 <= s
    @assert 0 <= b
    return z*lerch(-z,s,1.0,b,rtol=rtol)
end 

end