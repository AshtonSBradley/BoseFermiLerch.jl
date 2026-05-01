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
    budget_us::Float64
end

const SPEED_CASES = [
    SpeedCase("complete_well_damped", "lerch", 0.5, 2.0, 1.0, 0.0, 400.0),
    SpeedCase("complete_near_one", "lerch", 0.999999, 2.0, 1.0, 0.0, 1500.0),
    SpeedCase("incomplete_small_b", "lerch", 0.5, 2.0, 1.0, 0.05, 1500.0),
    SpeedCase("incomplete_near_one_small_b", "lerch", 0.999, 2.0, 1.0, 0.05, 1500.0),
    SpeedCase("incomplete_just_above_one_small_b", "lerch", 1.01, 2.0, 1.0, 0.05, 1500.0),
    SpeedCase("incomplete_large_b", "lerch", 0.5, 2.0, 1.0, 1.0, 1500.0),
    SpeedCase("real_cut_window", "lerch", 1.5, 2.0, 1.0, 1.0, 1500.0),
    SpeedCase("fermi_complete_equiv", "lerch", -0.8, 2.0, 1.0, 0.0, 600.0),
    SpeedCase("complex_complete", "lerch", 0.6 + 0.2im, 2.0, 1.0, 0.0, 2500.0),
    SpeedCase("complex_incomplete", "lerch", 0.6 + 0.2im, 2.0, 1.0, 0.5, 3000.0),
    SpeedCase("large_order_complex", "lerch", -5.0 - 2.1im, 100.0, 1.0, 0.0, 100.0),
    SpeedCase("very_large_order_complex", "lerch", -5.0 - 2.1im, 2000.0, 1.0, 0.0, 20.0),
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

    median_us = median(times) * 1e6
    minimum_us = minimum(times) * 1e6
    memory_bytes = round(Int, median(byte_counts))
    allocs = round(Int, median(alloc_counts))
    value = lerch(z, s, a, b)

    return (; median_us, minimum_us, memory_bytes, allocs, value)
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
    println(rpad("case", 38), rpad("args (z,s,a,b)", 34), lpad("median us", 12), lpad("min us", 12), lpad("budget", 10), lpad("allocs", 10), lpad("bytes", 12), "  status")
    println("-"^138)
end

function main()
    print_header()
    failures = String[]

    for sc in SPEED_CASES
        result = run_case(sc)
        args = "(" * join((format_arg(sc.z), format_arg(sc.s), format_arg(sc.a), format_arg(sc.b)), ", ") * ")"
        status = result.median_us <= sc.budget_us ? "ok" : "over"

        println(
            rpad(sc.name, 38),
            rpad(args, 34),
            lpad(@sprintf("%.3f", result.median_us), 12),
            lpad(@sprintf("%.3f", result.minimum_us), 12),
            lpad(@sprintf("%.3f", sc.budget_us), 10),
            lpad(string(result.allocs), 10),
            lpad(string(result.memory_bytes), 12),
            "  ",
            status,
        )

        if ENFORCE && status != "ok"
            push!(failures, "$(sc.name): median $(round(result.median_us; digits = 3)) us > budget $(round(sc.budget_us; digits = 3)) us")
        end
    end

    if ENFORCE && !isempty(failures)
        error("speed regression failures:\n" * join(failures, "\n"))
    end
end

main()
