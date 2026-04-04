using Printf

const _pkg_root = normpath(joinpath(@__DIR__, ".."))
if _pkg_root ∉ LOAD_PATH
    pushfirst!(LOAD_PATH, _pkg_root)
end

using BenchmarkTools
using BoseFermiLerch

const BENCH_SECONDS = 0.08
const WARMUP_SECONDS = 0.02

const REAL_COMPLETE_Z = [0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999]
const ORDERS = [1.5, 2.0, 3.0, 5.0]
const SHIFTS = [1.0, 2.0, 5.0]
const INCOMPLETE_B = [1e-3, 0.05, 0.5, 1.0]
const REAL_CUT_CASES = [
    (1.001, 0.01),
    (1.01, 0.05),
    (1.05, 0.1),
    (1.5, 1.0),
]
const COMPLEX_CASES = [
    (0.6 + 0.2im, 0.0),
    (0.95 + 0.05im, 0.0),
    (0.6 + 0.2im, 0.05),
    (0.95 + 0.05im, 0.5),
]

function bench_expr(f)
    trial = @benchmarkable $f seconds = BENCH_SECONDS
    tune!(trial)
    run(trial; seconds = WARMUP_SECONDS)
    return run(trial; seconds = BENCH_SECONDS)
end

function bench_lerch(z, s, a, b; rtol = 1e-9)
    lerch(z, s, a, b; rtol = rtol)
    trial = @benchmarkable lerch($z, $s, $a, $b; rtol = $rtol) seconds = BENCH_SECONDS
    tune!(trial)
    run(trial; seconds = WARMUP_SECONDS)
    result = run(trial; seconds = BENCH_SECONDS)
    return median(result).time / 1e6
end

function bench_complete_backends(z, s, a; rtol = 1e-10)
    vs = BoseFermiLerch.lerch_series(z, s, a, 0; rtol = rtol)
    vc = BoseFermiLerch.lerch_complete_contour(z, s, a; rtol = rtol)

    series_trial = @benchmarkable BoseFermiLerch.lerch_series($z, $s, $a, 0; rtol = $rtol) seconds = BENCH_SECONDS
    contour_trial = @benchmarkable BoseFermiLerch.lerch_complete_contour($z, $s, $a; rtol = $rtol) seconds = BENCH_SECONDS
    tune!(series_trial)
    tune!(contour_trial)
    run(series_trial; seconds = WARMUP_SECONDS)
    run(contour_trial; seconds = WARMUP_SECONDS)
    series_result = run(series_trial; seconds = BENCH_SECONDS)
    contour_result = run(contour_trial; seconds = BENCH_SECONDS)

    return (
        series_ms = median(series_result).time / 1e6,
        contour_ms = median(contour_result).time / 1e6,
        rel_diff = abs(vs - vc) / abs(vc),
    )
end

function print_section(title)
    println()
    println(title)
    println("="^length(title))
end

function main()
    println("BoseFermiLerch dense regime matrix")
    println("benchmark seconds per sample: ", BENCH_SECONDS)

    print_section("Complete real grid (b = 0)")
    println(rpad("z", 10), rpad("s", 8), rpad("a", 8), lpad("lerch ms", 12))
    println("-"^38)
    for z in REAL_COMPLETE_Z, s in ORDERS, a in SHIFTS
        t = bench_lerch(z, s, a, 0.0)
        println(rpad(@sprintf("%.5f", z), 10), rpad(@sprintf("%.3f", s), 8), rpad(@sprintf("%.1f", a), 8), lpad(@sprintf("%.6f", t), 12))
    end

    print_section("Incomplete real grid (0 < z < 1, b > 0)")
    println(rpad("z", 10), rpad("s", 8), rpad("a", 8), rpad("b", 10), lpad("lerch ms", 12))
    println("-"^48)
    for z in (0.5, 0.99, 0.999), s in ORDERS, a in SHIFTS, b in INCOMPLETE_B
        t = bench_lerch(z, s, a, b)
        println(rpad(@sprintf("%.5f", z), 10), rpad(@sprintf("%.3f", s), 8), rpad(@sprintf("%.1f", a), 8), rpad(@sprintf("%.5f", b), 10), lpad(@sprintf("%.6f", t), 12))
    end

    print_section("Real cut window (1 < z < exp(b))")
    println(rpad("z", 10), rpad("s", 8), rpad("a", 8), rpad("b", 10), lpad("lerch ms", 12))
    println("-"^48)
    for (z, b) in REAL_CUT_CASES, s in (1.5, 2.0, 3.0), a in (1.0, 2.0)
        t = bench_lerch(z, s, a, b)
        println(rpad(@sprintf("%.5f", z), 10), rpad(@sprintf("%.3f", s), 8), rpad(@sprintf("%.1f", a), 8), rpad(@sprintf("%.5f", b), 10), lpad(@sprintf("%.6f", t), 12))
    end

    print_section("Complex representative cases")
    println(rpad("z", 24), rpad("s", 8), rpad("a", 8), rpad("b", 10), lpad("lerch ms", 12))
    println("-"^62)
    for (z, b) in COMPLEX_CASES, s in (1.5, 2.0, 3.0), a in (1.0, 2.0)
        t = bench_lerch(z, s, a, b)
        println(rpad(string(z), 24), rpad(@sprintf("%.3f", s), 8), rpad(@sprintf("%.1f", a), 8), rpad(@sprintf("%.5f", b), 10), lpad(@sprintf("%.6f", t), 12))
    end

    print_section("Complete backend crossover (series vs contour)")
    println(rpad("z", 10), rpad("s", 8), rpad("a", 8), lpad("series ms", 12), lpad("contour ms", 13), lpad("rel diff", 12))
    println("-"^63)
    for z in REAL_COMPLETE_Z, s in (1.5, 2.0, 3.0), a in (1.0,)
        result = try
            bench_complete_backends(z, s, a)
        catch err
            println(rpad(@sprintf("%.5f", z), 10), rpad(@sprintf("%.3f", s), 8), rpad(@sprintf("%.1f", a), 8), lpad("series fail", 12), lpad("-", 13), lpad("-", 12))
            continue
        end
        println(
            rpad(@sprintf("%.5f", z), 10),
            rpad(@sprintf("%.3f", s), 8),
            rpad(@sprintf("%.1f", a), 8),
            lpad(@sprintf("%.6f", result.series_ms), 12),
            lpad(@sprintf("%.6f", result.contour_ms), 13),
            lpad(@sprintf("%.3e", result.rel_diff), 12),
        )
    end
end

main()
