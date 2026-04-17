using Printf
using Statistics

const _pkg_root = normpath(joinpath(@__DIR__, ".."))
if _pkg_root ∉ LOAD_PATH
    pushfirst!(LOAD_PATH, _pkg_root)
end

using Chairmarks
using BoseFermiLerch

struct SpeedCase
    name::String
    family::String
    z
    s
    a
    b
    budget_ms::Float64
end

const SPEED_CASES = [
    SpeedCase("complete_well_damped", "lerch", 0.5, 2.0, 1.0, 0.0, 0.40),
    SpeedCase("complete_near_one", "lerch", 0.999999, 2.0, 1.0, 0.0, 1.50),
    SpeedCase("incomplete_small_b", "lerch", 0.5, 2.0, 1.0, 0.05, 1.50),
    SpeedCase("incomplete_large_b", "lerch", 0.5, 2.0, 1.0, 1.0, 1.50),
    SpeedCase("real_cut_window", "lerch", 1.5, 2.0, 1.0, 1.0, 1.50),
    SpeedCase("fermi_complete_equiv", "lerch", -0.8, 2.0, 1.0, 0.0, 0.60),
    SpeedCase("complex_complete", "lerch", 0.6 + 0.2im, 2.0, 1.0, 0.0, 2.50),
    SpeedCase("complex_incomplete", "lerch", 0.6 + 0.2im, 2.0, 1.0, 0.5, 3.00),
]

const WARMUP_SECONDS = 0.05
const BENCH_SECONDS = 0.20
const ENFORCE = get(ENV, "BFL_ENFORCE_SPEED", "0") == "1"

function run_case(sc::SpeedCase)
    z, s, a, b = sc.z, sc.s, sc.a, sc.b

    # Warm up the compilation/runtime path before timing.
    lerch(z, s, a, b)

    run = @be lerch(z, s, a, b) seconds = BENCH_SECONDS
    times = getproperty.(run.samples, :time)
    alloc_counts = getproperty.(run.samples, :allocs)
    byte_counts = getproperty.(run.samples, :bytes)

    median_ms = median(times) * 1e3
    minimum_ms = minimum(times) * 1e3
    memory_bytes = round(Int, median(byte_counts))
    allocs = round(Int, median(alloc_counts))
    value = lerch(z, s, a, b)

    return (; median_ms, minimum_ms, memory_bytes, allocs, value)
end

function format_arg(x)
    if x isa Complex
        return string(x)
    elseif x isa AbstractFloat
        return @sprintf("%.6g", x)
    else
        return string(x)
    end
end

function print_header()
    println("BoseFermiLerch speed regression report")
    println("mode: ", ENFORCE ? "strict budgets" : "report only")
    println("benchmark seconds per case: ", BENCH_SECONDS)
    println()
    println(rpad("case", 24), rpad("args (z,s,a,b)", 34), lpad("median ms", 12), lpad("min ms", 12), lpad("budget", 10), lpad("allocs", 10), lpad("bytes", 12), "  status")
    println("-"^126)
end

function main()
    print_header()
    failures = String[]

    for sc in SPEED_CASES
        result = run_case(sc)
        args = "(" * join((format_arg(sc.z), format_arg(sc.s), format_arg(sc.a), format_arg(sc.b)), ", ") * ")"
        status = result.median_ms <= sc.budget_ms ? "ok" : "over"

        println(
            rpad(sc.name, 24),
            rpad(args, 34),
            lpad(@sprintf("%.3f", result.median_ms), 12),
            lpad(@sprintf("%.3f", result.minimum_ms), 12),
            lpad(@sprintf("%.3f", sc.budget_ms), 10),
            lpad(string(result.allocs), 10),
            lpad(string(result.memory_bytes), 12),
            "  ",
            status,
        )

        if ENFORCE && status != "ok"
            push!(failures, "$(sc.name): median $(round(result.median_ms; digits = 3)) ms > budget $(round(sc.budget_ms; digits = 3)) ms")
        end
    end

    if ENFORCE && !isempty(failures)
        error("speed regression failures:\n" * join(failures, "\n"))
    end
end

main()
