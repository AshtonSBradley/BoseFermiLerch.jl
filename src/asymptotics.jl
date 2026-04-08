const Z1_REAL_ASYMP_MU_CUTOFF = 0.05
const Z1_REAL_ASYMP_MINTERMS = 6
const Z1_REAL_ASYMP_MAXTERMS = 64

"""
    _near_integer_order(s::T) where {T<:Real}

Return `(is_integer_like, n)` for a positive real order `s`, using a tolerance
scaled to the working precision so the near-integer branch is handled
deliberately rather than by accidental `Inf - Inf` cancellation.
"""
function _near_integer_order(s::T) where {T<:Real}
    tol = T(32) * eps(float(s))
    n = round(Int, s)
    return abs(s - T(n)) <= tol, n
end

"""
    _harmonic_number(::Type{T}, n) where {T<:Real}

Compute `H_n = sum_{k=1}^n 1/k` in the target arithmetic type.
"""
function _harmonic_number(::Type{T}, n::Integer) where {T<:Real}
    total = zero(T)
    for k in 1:n
        total += inv(T(k))
    end
    return total
end

"""
    _should_use_complete_z1_real_asymptotic(z, s, a)

Conservative dispatch predicate for the complete positive-real near-degenerate
Lerch/Bose regime. The patch is intentionally limited to the complete case with
`a = 1`, real `0 < z < 1`, real `s > 0`, and `|log(z)| <= 0.05`.
"""
function _should_use_complete_z1_real_asymptotic(z, s, a)
    a == one(a) || return false
    isreal(z) && isreal(s) || return false
    zr = real(z)
    sr = real(s)
    zero(zr) < zr < one(zr) || return false
    sr > zero(sr) || return false
    return abs(log(zr)) <= Z1_REAL_ASYMP_MU_CUTOFF
end

"""
    _polylog_z1_regular_sum_noninteger(mu::T, s::T; atol, rtol, maxterms, minterms)
    where {T<:Real}

Accumulate the regular zeta-series part of the Robinson expansion for
`Li_s(exp(mu))` when `s` is not near a positive integer:

`sum_{k >= 0} zeta(s-k) * mu^k / k!`.

The series is asymptotic, so the implementation enforces a minimum number of
terms, stops when the terms are small relative to the partial sum, and falls
back to the smallest-term principle if the tail starts growing.
"""
function _polylog_z1_regular_sum_noninteger(
    mu::T,
    s::T;
    atol::T,
    rtol::T,
    maxterms::Int,
    minterms::Int,
) where {T<:Real}
    mu_pow = one(T)
    inv_fact = one(T)
    total = zeta(s)
    best_value = total
    best_term_abs = abs(total)
    prev_term_abs = best_term_abs
    grow_streak = 0

    for k in 1:maxterms
        mu_pow *= mu
        inv_fact /= T(k)
        term = zeta(s - T(k)) * mu_pow * inv_fact
        total += term

        term_abs = abs(term)
        if term_abs <= best_term_abs
            best_term_abs = term_abs
            best_value = total
        end

        if k >= minterms
            tol = max(atol, rtol * max(abs(total), one(T)))
            if term_abs <= tol
                return (; value = total, nterms = k + 1, stop_reason = :relative_tolerance)
            end

            if term_abs > prev_term_abs
                grow_streak += 1
                if grow_streak >= 2
                    return (; value = best_value, nterms = k + 1, stop_reason = :smallest_term)
                end
            else
                grow_streak = 0
            end
        end

        prev_term_abs = term_abs
    end

    return (; value = best_value, nterms = maxterms + 1, stop_reason = :maxterms)
end

"""
    _polylog_z1_integer_expansion(mu::T, n::Int; atol, rtol, maxterms, minterms)
    where {T<:Real}

Evaluate the logarithmic continuation for `Li_n(exp(mu))` with positive integer
`n` and real `mu < 0`:

`Li_n(e^mu) = sum_{k>=0, k!=n-1} zeta(n-k) mu^k/k!
            + mu^(n-1)/(n-1)! * (H_{n-1} - log(-mu))`

This avoids the invalid `Inf - Inf` cancellation between the
`Gamma(1-s) * (-mu)^(s-1)` singular term and the `zeta(1)` term of the regular
series.
"""
function _polylog_z1_integer_expansion(
    mu::T,
    n::Int;
    atol::T,
    rtol::T,
    maxterms::Int,
    minterms::Int,
) where {T<:Real}
    n == 1 && return (; value = -log1p(-exp(mu)), nterms = 0, stop_reason = :exact_log)

    inv_fact_log = one(T)
    for k in 1:(n - 1)
        inv_fact_log /= T(k)
    end
    log_term = mu^(n - 1) * inv_fact_log * (_harmonic_number(T, n - 1) - log(-mu))

    total = log_term + zeta(T(n))
    best_value = total
    best_term_abs = abs(zeta(T(n)))
    prev_term_abs = best_term_abs
    grow_streak = 0
    mu_pow = one(T)
    inv_fact = one(T)

    for k in 1:maxterms
        mu_pow *= mu
        inv_fact /= T(k)

        if k == n - 1
            continue
        end

        term = zeta(T(n - k)) * mu_pow * inv_fact
        total += term

        term_abs = abs(term)
        if term_abs <= best_term_abs
            best_term_abs = term_abs
            best_value = total
        end

        if k >= minterms
            tol = max(atol, rtol * max(abs(total), one(T)))
            if term_abs <= tol
                return (; value = total, nterms = k + 1, stop_reason = :relative_tolerance)
            end

            if term_abs > prev_term_abs
                grow_streak += 1
                if grow_streak >= 2
                    return (; value = best_value, nterms = k + 1, stop_reason = :smallest_term)
                end
            else
                grow_streak = 0
            end
        end

        prev_term_abs = term_abs
    end

    return (; value = best_value, nterms = maxterms + 1, stop_reason = :maxterms)
end

"""
    _bose_asymp_z1_real_stats(z::T, s::T; atol, rtol, maxterms, minterms) where {T<:Real}

Return the Robinson-style asymptotic evaluation data for the complete
positive-real Bose/polylog case near `z = 1`, where `mu = log(z) < 0`:

`Li_s(e^mu) = Gamma(1-s) * (-mu)^(s-1) + sum_{k>=0} zeta(s-k) mu^k/k!`

for noninteger `s`, with the positive-integer branch replaced by the exact
logarithmic continuation.
"""
function _bose_asymp_z1_real_stats(
    z::T,
    s::T;
    atol::T = zero(T),
    rtol::T = sqrt(eps(T)),
    maxterms::Int = Z1_REAL_ASYMP_MAXTERMS,
    minterms::Int = Z1_REAL_ASYMP_MINTERMS,
) where {T<:Real}
    _should_use_complete_z1_real_asymptotic(z, s, one(T)) ||
        throw(DomainError((z, s), "near-z=1 asymptotic patch is only defined for real 0 < z < 1, real s > 0, and a = 1"))

    mu = log(z)
    is_integer_like, n = _near_integer_order(s)

    if is_integer_like
        return merge((mu = mu,), _polylog_z1_integer_expansion(mu, n; atol = atol, rtol = rtol, maxterms = maxterms, minterms = minterms))
    end

    singular = gamma(one(T) - s) * (-mu)^(s - one(T))
    regular = _polylog_z1_regular_sum_noninteger(mu, s; atol = atol, rtol = rtol, maxterms = maxterms, minterms = minterms)
    return (; value = singular + regular.value, mu = mu, nterms = regular.nterms, stop_reason = regular.stop_reason)
end

"""
    _bose_asymp_z1_real(z::T, s::T; atol, rtol, maxterms, minterms) where {T<:Real}

Evaluate `g_s(z) = Li_s(z)` in the real complete near-`z = 1` Bose regime.
"""
function _bose_asymp_z1_real(
    z::T,
    s::T;
    atol::T = zero(T),
    rtol::T = sqrt(eps(T)),
    maxterms::Int = Z1_REAL_ASYMP_MAXTERMS,
    minterms::Int = Z1_REAL_ASYMP_MINTERMS,
) where {T<:Real}
    return _bose_asymp_z1_real_stats(z, s; atol = atol, rtol = rtol, maxterms = maxterms, minterms = minterms).value
end

function _bose_asymp_z1_real(
    z::Real,
    s::Real;
    atol::Real = zero(promote_type(typeof(float(z)), typeof(float(s)))),
    rtol::Real = sqrt(eps(promote_type(typeof(float(z)), typeof(float(s))))),
    maxterms::Int = Z1_REAL_ASYMP_MAXTERMS,
    minterms::Int = Z1_REAL_ASYMP_MINTERMS,
)
    T = promote_type(typeof(float(z)), typeof(float(s)), typeof(float(atol)), typeof(float(rtol)))
    return _bose_asymp_z1_real(T(z), T(s); atol = T(atol), rtol = T(rtol), maxterms = maxterms, minterms = minterms)
end

"""
    _lerch_complete_asymp_z1_real(z::T, s::T; atol, rtol, maxterms, minterms) where {T<:Real}

Evaluate the complete Lerch case `Phi(z, s, 1) = Li_s(z) / z` in the same
positive-real near-`z = 1` regime handled by `_bose_asymp_z1_real`.
"""
function _lerch_complete_asymp_z1_real(
    z::T,
    s::T;
    atol::T = zero(T),
    rtol::T = sqrt(eps(T)),
    maxterms::Int = Z1_REAL_ASYMP_MAXTERMS,
    minterms::Int = Z1_REAL_ASYMP_MINTERMS,
) where {T<:Real}
    return _bose_asymp_z1_real(z, s; atol = atol, rtol = rtol, maxterms = maxterms, minterms = minterms) / z
end

function _lerch_complete_asymp_z1_real(
    z::Real,
    s::Real;
    atol::Real = zero(promote_type(typeof(float(z)), typeof(float(s)))),
    rtol::Real = sqrt(eps(promote_type(typeof(float(z)), typeof(float(s))))),
    maxterms::Int = Z1_REAL_ASYMP_MAXTERMS,
    minterms::Int = Z1_REAL_ASYMP_MINTERMS,
)
    T = promote_type(typeof(float(z)), typeof(float(s)), typeof(float(atol)), typeof(float(rtol)))
    return _lerch_complete_asymp_z1_real(T(z), T(s); atol = T(atol), rtol = T(rtol), maxterms = maxterms, minterms = minterms)
end
