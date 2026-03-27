using Test
using SpecialFunctions

F0(x) = log1p(exp(x))

for x in -10:10
    @test fermi(1, exp(x)) ≈ F0(x)
end

@test fermi(1, 1) ≈ log(2)
@test fermi(2, 1) ≈ pi^2 / 12
@test fermi(2, 0.5, 0.3) ≈ 0.5 * lerch(-0.5, 2, 1, 0.3)
@test fermi(2, 0, 0.3) == 0
@test_throws DomainError fermi(0, 1)
@test_throws DomainError fermi(2, -2.0, 0.5)
