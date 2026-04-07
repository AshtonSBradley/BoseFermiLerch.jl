using Test
using SpecialFunctions

@test bose(1, 3) == bose(1, 3, 0)
@test bose(1, 5) == zeta(5)
@test bose(0.0, 1) == 0.0
@test bose(0.5, 1) ≈ log(2)
@test bose(0.999999999, 1) ≈ -log1p(-0.999999999)
@test bose(1 - 1e-12, 1) ≈ -log(1e-12)
@test bose(0.3 + 0.2im, 1) ≈ -log1p(-(0.3 + 0.2im))
@test bose(0.5, 2, 0.3) ≈ 0.5 * lerch(0.5, 2, 1, 0.3)
@test bose(0.999999, 2, 0.0) ≈ 0.999999 * lerch(0.999999, 2, 1, 0.0)
@test bose(0, 2, 0.3) == 0
@test_throws DomainError bose(0.5, 0)
@test_throws DomainError bose(2.0, 2, 0.5)
@test_throws ArgumentError bose(0.5, 2; rtol = 0)
@test_throws ArgumentError bose(0.5, 2; rtol = -1e-8)

#TODO: tidy these
# ## Lerchphi test 
# z,s,a,b =.5,1,1.,0.

# ## bose
# lerch(z,s,a,b)
# @btime lerch(z,s,a,b)

# bose(z,s,b)
# z*lerch(z,s,a,b)

# @btime bose(z,s,b)
# @btime z*lerch(z,s,a,b)

# # fermi
# @btime fermi(z,s,b)
# @btime z*lerch(-z,s,a,b)

# z,s,a,b =0.5,1.0,1.,0.
# z*lerch(z,s,a,b)
# bose(z,s,b)

# z,s,a,b = 1.0,2.0,1.,0.
# @btime z*lerch(z,s,a,b)
# @btime bose(z,s,b)

## simple timing

    # @btime bose(1,3/2,0.)
    # @btime zeta(3/2)
    # @btime bose(.5,3,.5)
    # exp(.5)
    # zeta(3)
    # bose(.99,0.1,.1)

    # # note only getting 5 digits for log(2)
    # @test bose(1,3)==bose(1,3,0)
    # @test bose(1,5)==zeta(5)
    # @test bose(.5,1) ≈ log(2)
