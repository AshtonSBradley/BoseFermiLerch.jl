
@test bose(3,1) == bose(3,1,0)
@test bose(5,1) == zeta(5)
@test bose(1,.5) ≈ log(2)

#TODO: tidy these
# ## Lerchphi test 
# z,s,a,b =.5,1,1.,0.

# ## bose
# lerch(z,s,a,b)
# @btime lerch(z,s,a,b)

# bose(s,z,b)
# z*lerch(z,s,a,b)

# @btime bose(s,z,b)
# @btime z*lerch(z,s,a,b)

# # fermi
# @btime fermi(s,z,b)
# @btime z*lerch(-z,s,a,b)

# z,s,a,b =0.5,1.0,1.,0.
# z*lerch(z,s,a,b)
# bose(s,z,b)

# z,s,a,b = 1.0,2.0,1.,0.
# @btime z*lerch(z,s,a,b)
# @btime bose(s,z,b)

## simple timing

    # @btime bose(3/2,1,0.)
    # @btime zeta(3/2)
    # @btime bose(3,.5,.5)
    # exp(.5)
    # zeta(3)
    # bose(0.1,.99,.1)

    # # note only getting 5 digits for log(2)
    # @test bose(3,1)==bose(3,1,0)
    # @test bose(5,1)==zeta(5)
    # @test bose(1,.5) ≈ log(2)