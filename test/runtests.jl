using BoseFermi
using SpecialFunctions
using Test

≈(x,y)=isapprox(x,y,rtol=1e-4)
@testset "bose tests" begin include("bose_tests.jl") end

