
@test bose(3,1) == bose(3,1,0)
@test bose(5,1) == zeta(5)
@test bose(1,.5) ≈ log(2)
