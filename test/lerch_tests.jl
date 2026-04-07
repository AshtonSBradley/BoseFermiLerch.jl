using Test
using SpecialFunctions
using QuadGK

lerch_neg_real_ref(z, s, a; rtol = 1e-8) =
    quadgk(t -> exp(-a * t) * t^(s - 1) / ((1 - z * exp(-t)) * gamma(s)), 0.0, Inf; rtol = rtol)[1]

@test lerch(1.0, 2, 1, 0, rtol = 1e-9) ≈ zeta(2)
@test lerch(0.5, 1, 1, 0, rtol = 1e-10) ≈ -log(1 - 0.5) / 0.5
@test lerch(-0.5, 2, 1, 0, rtol = 1e-10) ≈ sum((-0.5)^n / (n + 1)^2 for n in 0:200)
@test lerch(-10.0, 3 / 2, 1, 0; rtol = 1e-7) ≈ lerch_neg_real_ref(-10.0, 3 / 2, 1; rtol = 1e-7)
@test lerch(0.5, 2, 1, 0.3, rtol = 1e-10) ≈ sum(0.5^n * gamma(2, 0.3 * (n + 1)) / (n + 1)^2 for n in 0:200)
@test lerch(0.8, 2, 1, 0.05, rtol = 1e-10) ≈ sum(0.8^n * gamma(2, 0.05 * (n + 1)) / (n + 1)^2 for n in 0:300)
@test lerch(1.0, 2, 1, 0.3, rtol = 1e-10) ≈ sum(gamma(2, 0.3 * (n + 1)) / (n + 1)^2 for n in 0:300)
@test isapprox(lerch(0.999999, 2, 1, 0, rtol = 1e-10), zeta(2) / 0.999999; atol = 1e-4)
@test lerch(0.999, 2, 1, 0.3, rtol = 1e-9) ≈ sum(0.999^n * gamma(2, 0.3 * (n + 1)) / (n + 1)^2 for n in 0:500)
@test lerch(1.5, 2, 1, 1.0, rtol = 1e-9) ≈ sum(1.5^n * gamma(2, n + 1) / (n + 1)^2 for n in 0:500)
let b = 0.1, z = exp(b) - 1e-6
    ref = quadgk(t -> exp(-t) * t / (1 - z * exp(-t)), b, Inf; rtol = 1e-9)[1]
    @test lerch(z, 2, 1, b; rtol = 1e-9) ≈ ref
end
let z = 0.7 + 0.2im, s = 2, a = 0.5
    @test lerch(z, s, a, 0; rtol = 1e-9) ≈ z * lerch(z, s, a + 1, 0; rtol = 1e-9) + a^(-s)
end
@test_throws DomainError lerch(1.5, 2, 1, 0)
@test_throws DomainError lerch(3.0, 2, 1, 0.5)
@test_throws DomainError lerch(0.5, 2, 0.0, 0.1)
@test_throws DomainError lerch(0.5, 2, -1.0, 0.1)
@test_throws DomainError lerch(0.5, 2, 1, -0.1)
@test_throws ArgumentError lerch(0.5, 2, 1, 0.1; rtol = 0)
