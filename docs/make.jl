
import Pkg; Pkg.add("Documenter")

using Documenter, BoseFermi

makedocs(sitename="BoseFermi.jl",
pages = [
"Overview" => "index.md",
"Bose" => "bose.md",
"Fermi" => "fermi.md",
"Lerch" => "lerch.md"
],
format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/AshtonSBradley/BoseFermi.jl.git",
)