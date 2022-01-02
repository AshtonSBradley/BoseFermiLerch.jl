using BoseFermi

@test bose(1,3)==bose(1,3,0)
@test bose(1,5,)==zeta(5)
@test bose(.5,1)≈log(2)
