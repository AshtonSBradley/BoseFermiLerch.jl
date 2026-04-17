using BoseFermiLerch
using SpecialFunctions
using Test

≈(x,y)=isapprox(x,y,rtol=1e-4)
@testset "lerch tests" begin include("lerch_tests.jl") end
@testset "bose tests" begin include("bose_tests.jl") end
@testset "fermi tests" begin include("fermi_tests.jl") end
@testset "bose near one tests" begin include("test_bose_near_one.jl") end
@testset "asymptotic tests" begin include("test_asymptotics.jl") end
