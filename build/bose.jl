gamma_inc_Q(p,x)=gamma_inc(p,x,0)[2]

"""
`bose(ν,z,y)`

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
This implementation requires the normalized incomplete gamma function provided by `GSL.jl`.
"""
function bose(ν,z,y=0.,rtol=1e-6,atol=rtol)
    @assert 0 <= ν
    @assert 0 <= y
    (z == 1 && y == 0) && return zeta(ν)
    (ν == 1 && z == 0.5 ) && return log(2)
    k = 1
    Sk = z^k/k^ν*gamma_inc_Q(ν,k*y)
    Skplus1 = Sk + z^(k+1)/(k+1)^ν*gamma_inc_Q(ν,(k+1)*y)
    while !isapprox(Sk,Skplus1,rtol=rtol,atol=atol)
        k += 1
        Sk = Skplus1
        Skplus1 += z^(k+1)/(k+1)^ν*gamma_inc_Q(ν,(k+1)*y)
    end
    return Sk
end