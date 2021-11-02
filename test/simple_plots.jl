using BoseFermi
using Plots, LaTeXStrings

## 
ϵ = 0.1
μ = LinRange(-.5,-.1,100)
z = exp.(μ)

##
ν = 2
plot(z,bose.(ν,z))
