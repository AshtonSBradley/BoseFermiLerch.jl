using Test
using SpecialFunctions
using QuadGK

F0(x) = log1p(exp(x))

for x in -10:10
    @test fermi(exp(x), 1) ≈ F0(x)
end

@test fermi(1, 1) ≈ log(2)
@test fermi(1, 2) ≈ pi^2 / 12
@test fermi(0.999999999, 1) ≈ log1p(0.999999999)
@test fermi(0.3 + 0.2im, 1) ≈ log1p(0.3 + 0.2im)
@test fermi(0.5, 2, 0.3) ≈ 0.5 * lerch(-0.5, 2, 1, 0.3)
@test fermi(-0.999999, 2, 0.0) ≈ -0.999999 * lerch(0.999999, 2, 1, 0.0)
@test fermi(0, 2, 0.3) == 0
@test_throws DomainError fermi(1, 0)
@test_throws DomainError fermi(-2.0, 2, 0.5)
@test_throws ArgumentError fermi(0.5, 2; rtol = 0)

fd_ref(s, z; rtol = 1e-8) = z * quadgk(t -> exp(-t) * t^(s - 1) / ((1 + z * exp(-t)) * gamma(s)), 0.0, Inf; rtol = rtol)[1]

@test fermi(exp(30), 3 / 2) ≈ fd_ref(3 / 2, exp(30))
@test fermi(exp(30), 5 / 2) ≈ fd_ref(5 / 2, exp(30))
@test fermi(exp(50), 3 / 2; rtol = 1e-7) ≈ fd_ref(3 / 2, exp(50); rtol = 1e-7)
@test fermi(exp(50), 5 / 2; rtol = 1e-7) ≈ fd_ref(5 / 2, exp(50); rtol = 1e-7)

let a = 3 / 2
    entropies = map([25, 30, 35, 40, 45, 50]) do x
        z = exp(x)
        f1 = fermi(z, a)
        f2 = fermi(z, a + 1)
        (a + 1) * f2 / f1 - log(z)
    end
    @test all(>(0), entropies)
    @test all(diff(entropies) .< 0)
end
