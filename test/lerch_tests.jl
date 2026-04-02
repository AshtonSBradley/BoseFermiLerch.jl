using Test
using SpecialFunctions

@test lerch(1.0, 2, 1, 0, rtol = 1e-9) ≈ zeta(2)
@test lerch(0.5, 1, 1, 0, rtol = 1e-10) ≈ -log(1 - 0.5) / 0.5
@test lerch(0.5, 2, 1, 0.3, rtol = 1e-10) ≈ sum(0.5^n * gamma(2, 0.3 * (n + 1)) / (n + 1)^2 for n in 0:200)
@test isapprox(lerch(0.999999, 2, 1, 0, rtol = 1e-10), zeta(2) / 0.999999; atol = 1e-4)
@test_throws DomainError lerch(3.0, 2, 1, 0.5)
@test_throws DomainError lerch(0.5, 2, 0.0, 0.1)
@test_throws DomainError lerch(0.5, 2, -1.0, 0.1)
