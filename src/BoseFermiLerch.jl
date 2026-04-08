module BoseFermiLerch

import SpecialFunctions: zeta, gamma, loggamma
using QuadGK

include("asymptotics.jl")

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

function check_rtol(rtol)
    rtol > zero(rtol) || throw(ArgumentError("rtol must be positive"))
end

function is_valid_incomplete_real_cut_point(z, b)
    return isreal(z) && one(real(z)) < real(z) < exp(b)
end

function series_cutoff(z, s)
    if isreal(z) && real(z) >= 0
        return isreal(s) && real(s) >= 2 ? 0.98 : 0.9
    end
    return 0.8
end

function should_use_series(z, s)
    abs(z) < 1 || return false
    return abs(z) <= series_cutoff(z, s)
end

function should_use_real_integral(z)
    return isreal(z) && real(z) < 0
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

function lerch_complete_contour(z, s, a; rtol = 1e-8)
    if isreal(z) && real(z) > 1
        throw(DomainError(z, "complete Lerch principal branch is undefined on the real interval (1, Inf)"))
    end

    if z == one(z)
        return zeta(s, a)
    end

    if real(a) < 2
        return z * lerch_complete_contour(z, s, a + 1; rtol = rtol) + a^(-s)
    end

    g(t) = t^(s - 1) * exp(-a * t) / (1 - z * exp(-t))
    h(t) = (-t)^(s - 1) * exp(-a * t) / (1 - z * exp(-t))
    L = log(complex(z))

    if isinteger(s) && real(s) >= 1
        if abs(imag(L)) < 0.25 && real(L) >= 0
            if imag(z) <= 0
                I = quadgk(g, 0.0, 1im, 1im + abs(L) + 1, abs(L) + 1, Inf; rtol = rtol)[1]
            else
                I = quadgk(g, 0.0, -1im, -1im + abs(L) + 1, abs(L) + 1, Inf; rtol = rtol)[1]
            end
        else
            I = quadgk(g, 0.0, Inf; rtol = rtol)[1]
        end
        return exp(-loggamma(s)) * I
    end

    if real(L) < -0.5
        residue = 0.0
        c = min(abs(real(L)) / 2, 1.0)
        left = right = top = c
    elseif abs(imag(L)) > 0.5
        residue = 0.0
        c = min(abs(imag(L)) / 2, 1.0)
        left = right = top = c
    else
        residue = (-L)^s / (L * z^a)
        left = max(0.0, -real(L)) + 1
        top = abs(imag(L)) + 1
        right = abs(L) + 1
    end

    is_real_case = isreal(z) && real(z) < 1 && isreal(s) && isreal(a) && real(a) > 0
    w = Complex(-1.0)^(s - 1)
    I = 0.0 + 0.0im

    if is_real_case
        I += 2im * imag(quadgk(g, right, right + top * im; rtol = rtol)[1] / w)
        I += 2im * imag(quadgk(g, right + top * im, -left + top * im; rtol = rtol)[1] / w)
        I += 2im * imag(quadgk(h, -left + top * im, -left; rtol = rtol)[1])
        I += quadgk(g, right, Inf; rtol = rtol)[1] * (w - inv(w))
    else
        I += quadgk(g, right, right + top * im; rtol = rtol)[1] / w
        I += quadgk(g, right + top * im, -left + top * im; rtol = rtol)[1] / w
        I += quadgk(h, -left + top * im, -left - top * im; rtol = rtol)[1]
        I += quadgk(g, -left - top * im, right - top * im; rtol = rtol)[1] * w
        I += quadgk(g, right - top * im, right; rtol = rtol)[1] * w
        I += quadgk(g, right, Inf; rtol = rtol)[1] * (w - inv(w))
    end

    I = -gamma(1 - s) * (I / (2pi * im) + residue)
    return is_real_case ? real(I) : I
end

function lerch_lower_correction(z, s, a, b; rtol = 1e-8)
    return quadgk(t -> lerch_int(t, z, s, a), zero(b), b; rtol = rtol)[1]
end

function lerch_complete_real_integral(z, s, a; rtol = 1e-8)
    return quadgk(t -> lerch_int(t, z, s, a), zero(real(a)), Inf; rtol = rtol)[1]
end

"""
    lerch(z, s, a, b; rtol=1e-8)

Evaluate the upper incomplete Lerch transcendent

```math
\\Phi(z,s,a,b)=\\frac{1}{\\Gamma(s)}\\int_b^\\infty \\frac{t^{s-1}e^{-at}}{1-ze^{-t}}\\,dt.
```

The implementation uses a hybrid backend: a gamma-series path where it is broadly
valid and benchmark-favoured, a dedicated asymptotic patch for the complete
real positive near-`z = 1` `a = 1` regime, a direct real-axis integral for the
complete negative-real regime that stabilizes Fermi-Dirac evaluations at large
fugacity, and the contour integral as a robust fallback. For `b > 0`, the
incomplete function is evaluated as the complete value minus the finite-interval
correction on `[0, b]`, except for the narrow real interval `1 < z < exp(b)`
where the incomplete integral is valid but the complete principal branch sits
on its cut, so the upper integral is evaluated directly.
"""
function lerch(z, s, a, b; rtol = 1e-8)
    check_order(s)
    check_lower_limit(b)
    check_shift(a)
    check_rtol(rtol)

    if b == 0
        if _should_use_complete_z1_real_asymptotic(z, s, a)
            zr = real(z)
            sr = real(s)
            return _lerch_complete_asymp_z1_real(zr, sr; rtol = rtol)
        end

        if should_use_series(z, s)
            return lerch_series(z, s, a, b; rtol = rtol)
        end

        # On the negative real interval, the direct integral is regular and avoids
        # the severe cancellation that can appear in the contour formula for large
        # negative z, which directly affects high-fugacity Fermi evaluations.
        if should_use_real_integral(z)
            return lerch_complete_real_integral(z, s, a; rtol = rtol)
        end
        return lerch_complete_contour(z, s, a; rtol = rtol)
    end

    check_branch_cut(z, b)

    # At z = 1 the complete-minus-lower decomposition inherits the complete-case
    # endpoint singularity at t = 0, while the incomplete series is exponentially
    # convergent for b > 0.
    if z == one(z)
        return lerch_series(z, s, a, b; rtol = rtol)
    end

    if should_use_series(z, s)
        return lerch_series(z, s, a, b; rtol = rtol)
    end

    # For real 1 < z < exp(b), the incomplete integral is regular because the pole
    # sits below the lower limit, but the complete principal branch is on its cut.
    if is_valid_incomplete_real_cut_point(z, b)
        return quadgk(t -> lerch_int(t, z, s, a), b, Inf; rtol = rtol)[1]
    end

    return lerch_complete_contour(z, s, a; rtol = rtol) - lerch_lower_correction(z, s, a, b; rtol = rtol)
end

"""
`bose(z, ν, y)`

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
function bose(z, s, b = 0; rtol = 1e-8)
    check_order(s)
    check_lower_limit(b)
    check_rtol(rtol)
    if z == one(z) && b == 0
        return zeta(s)
    elseif s == one(s) && b == zero(b)
        return -log1p(-z)
    elseif b == zero(b) && _should_use_complete_z1_real_asymptotic(z, s, one(s))
        return _bose_asymp_z1_real(real(z), real(s); rtol = rtol)
    else
        return z * lerch(z, s, 1.0, b; rtol = rtol)
    end
end

"""
`fermi(z, ν, y)`

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
function fermi(z, s, b = 0; rtol = 1e-8)
    check_order(s)
    check_lower_limit(b)
    check_rtol(rtol)
    if s == one(s) && b == zero(b)
        return log1p(z)
    end
    return z * lerch(-z, s, 1.0, b; rtol = rtol)
end

end
