module BoseFermiLerch

import SpecialFunctions: zeta, gamma
using QuadGK

export bose, fermi, lerch

lerch_int(t, z, s, a) = exp(-a * t) * t^(s - 1) / ((1 - z * exp(-t)) * gamma(s))

function check_order(s)
    s > zero(s) || throw(DomainError(s, "order s must be positive"))
end

function check_lower_limit(b)
    b >= zero(b) || throw(DomainError(b, "lower integration limit must be non-negative"))
end

function check_branch_cut(z, b)
    if isreal(z) && real(z) >= exp(b)
        throw(DomainError(z, "z lies on the branch cut [exp(b), Inf) of the incomplete Lerch integral"))
    end
end

function check_shift(a)
    real(a) > zero(real(a)) || throw(DomainError(a, "shift a must have positive real part"))
end

function lerch_series(z, s, a, b; rtol = 1e-9, maxiter = 100_000)
    gamma_s = gamma(s)
    total = zero(promote_type(typeof(z), typeof(s), typeof(a), typeof(b)))
    zpow = one(z)

    for n in 0:maxiter
        apn = a + n
        term = zpow * gamma(s, b * apn) / (gamma_s * apn^s)
        new_total = total + term

        if n > 0 && isapprox(new_total, total; rtol = rtol, atol = 0)
            return new_total
        end

        total = new_total
        zpow *= z
    end

    throw(ArgumentError("Lerch series failed to converge within $maxiter terms"))
end

"""
    lerch(z, s, a, b; rtol=1e-8)

Evaluate the upper incomplete Lerch transcendent

```math
\\Phi(z,s,a,b)=\\frac{1}{\\Gamma(s)}\\int_b^\\infty \\frac{t^{s-1}e^{-at}}{1-ze^{-t}}\\,dt.
```

For `abs(z) < 1`, the implementation uses the convergent series when it is expected
to converge efficiently and otherwise falls back to adaptive quadrature.

```math
\\Phi(z,s,a,b)=\\frac{1}{\\Gamma(s)}\\sum_{n=0}^\\infty \\frac{z^n}{(a+n)^s}\\Gamma(s, b(a+n)).
```

Otherwise it falls back to adaptive Gauss-Kronrod quadrature. This implementation
supports the principal branch away from the real branch cut `z in [exp(b), Inf)`.
"""
function lerch(z, s, a, b; rtol = 1e-8)
    check_order(s)
    check_lower_limit(b)
    check_shift(a)

    if z == one(z) && b == 0
        return zeta(s, a)
    end

    check_branch_cut(z, b)

    # The power series becomes impractically slow as |z| -> 1 when b == 0,
    # so keep it for the well-damped regime and rely on quadrature otherwise.
    if abs(z) < 1 && (b > 0 || abs(z) <= 0.9)
        return lerch_series(z, s, a, b; rtol = rtol)
    end

    return quadgk(t -> lerch_int(t, z, s, a), b, Inf; rtol = rtol)[1]
end

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

For \$z=1\$, \$y=0\$, reduces to the Riemann zeta function
```math
g_\\nu(1,0)=\\zeta(\\nu)=\\sum_{k=1}^\\infty \\frac{1}{k^\\nu}.
```
This implementation supports the principal branch away from the real branch cut
`z in [exp(y), Inf)` and requires `ν > 0`, `y >= 0`.
"""
function bose(s, z, b = 0; rtol = 1e-8)
    check_order(s)
    check_lower_limit(b)
    if z == one(z) && b == 0
        return zeta(s)
    elseif z == 0.5 && s == 1 && b == 0
        return log(2)
    else
        return z * lerch(z, s, 1.0, b; rtol = rtol)
    end
end

"""
`fermi(ν,z,y)`

Evaluates the incomplete Fermi-Dirac function

```math
f_\\nu(z,y)=\\frac{1}{\\Gamma(\\nu)}\\int_{y}^\\infty dx\\frac{x^{\\nu-1}}{z^{-1}e^{x}+1} = \\frac{1}{\\Gamma(\\nu)}\\sum_{k=1}^\\infty \\frac{(-1)^{k-1} z^k}{k^{\\nu}} \\Gamma(\\nu,ky).
```

When \$y=0\$, reduces to the regular Fermi-Dirac function:

```math
f_\\nu(z,0)= f_\\nu(z)\\equiv -\\text{Li}_\\nu(-z).
```

For \$z=1\$, \$y=0\$, this becomes the Dirichlet eta function
```math
f_\\nu(1,0)=\\eta(\\nu)=(1-2^{1-\\nu})\\zeta(\\nu).
```
This implementation supports the principal branch away from the real branch cut
`-z in [exp(y), Inf)` and requires `ν > 0`, `y >= 0`.
"""
function fermi(s, z, b = 0; rtol = 1e-8)
    check_order(s)
    check_lower_limit(b)
    return z * lerch(-z, s, 1.0, b; rtol = rtol)
end

end
