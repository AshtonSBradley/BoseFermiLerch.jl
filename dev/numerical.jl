using QuadGK, BenchmarkTools
using SpecialFunctions

Φint(t,z,s,a) = 1/gamma(s)t^(s-1)*exp(-a*t)/(1-z*exp(-t))

Φ(z,s,a,b) = quadgk(t->Φint(t,z,s,a),b,Inf)

@btime Φ(0.5,1.0,1.,0.)

p = 0.5Φ(0.5,1.0,1.,0.)[1]