using Printf

const _pkg_root = normpath(joinpath(@__DIR__, ".."))
if _pkg_root ∉ LOAD_PATH
    pushfirst!(LOAD_PATH, _pkg_root)
end

using BenchmarkTools
using QuadGK
using SpecialFunctions
using BoseFermiLerch

const BENCH_SECONDS = 0.12
const WARMUP_SECONDS = 0.03
const ASYMP_Z = [0.99, 0.999, 0.9999, 0.99999]
const ASYMP_S = [0.5, 1.5, 2.5]

function big_bose_ref(z, s; rtol = big"1e-30")
    setprecision(256) do
        zb = BigFloat(z)
        sb = BigFloat(s)
        gb = gamma(sb)
        integrand(t) = exp(-t) * t^(sb - 1) * zb / ((1 - zb * exp(-t)) * gb)
        quadgk(integrand, big"0.0", big"1.0", big"4.0", big"16.0", big"64.0", Inf; rtol = rtol)[1]
    end
end

function bench_call(f::F) where {F<:Function}
    trial = @benchmarkable ($f)() seconds = BENCH_SECONDS
    tune!(trial)
    run(trial; seconds = WARMUP_SECONDS)
    result = run(trial; seconds = BENCH_SECONDS)
    return median(result).time / 1e6
end

println("Near-z=1 Bose asymptotic benchmark")
println("benchmark seconds per sample: ", BENCH_SECONDS)
println()
println(rpad("s", 8), rpad("z", 12), lpad("public ms", 12), lpad("asymp ms", 12), lpad("contour ms", 14), lpad("public relerr", 16), lpad("asymp relerr", 16), lpad("speedup", 12))
println("-"^102)

for s in ASYMP_S, z in ASYMP_Z
    public_value = bose(z, s; rtol = 1e-12)
    asymp_value = BoseFermiLerch._bose_asymp_z1_real(z, s; rtol = 1e-12)
    contour_value = z * BoseFermiLerch.lerch_complete_contour(z, s, 1.0; rtol = 1e-12)
    ref = Float64(big_bose_ref(z, s))

    public_ms = bench_call(() -> bose(z, s; rtol = 1e-12))
    asymp_ms = bench_call(() -> BoseFermiLerch._bose_asymp_z1_real(z, s; rtol = 1e-12))
    contour_ms = bench_call(() -> z * BoseFermiLerch.lerch_complete_contour(z, s, 1.0; rtol = 1e-12))

    public_relerr = abs(public_value - ref) / abs(ref)
    asymp_relerr = abs(asymp_value - ref) / abs(ref)
    _ = contour_value

    println(
        rpad(@sprintf("%.3f", s), 8),
        rpad(@sprintf("%.5f", z), 12),
        lpad(@sprintf("%.6f", public_ms), 12),
        lpad(@sprintf("%.6f", asymp_ms), 12),
        lpad(@sprintf("%.6f", contour_ms), 14),
        lpad(@sprintf("%.3e", public_relerr), 16),
        lpad(@sprintf("%.3e", asymp_relerr), 16),
        lpad(@sprintf("%.2fx", contour_ms / asymp_ms), 12),
    )
end
