using Test
F0(x) = log(1+exp(x))

for x in -10:10
    @test fermi(1,exp(x)) ≈ F0(x)
end

