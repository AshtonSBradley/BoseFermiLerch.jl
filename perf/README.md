# Performance Regression Harness

This folder contains a small benchmark matrix for `lerch(z, s, a, b)` that probes:

- well-damped complete values
- complete values near `z = 1`, which now exercise the near-one asymptotic patch
- incomplete values with small and large `b`
- the real window `1 < z < exp(b)` where the incomplete integral is valid
- complex complete and incomplete inputs

These timings are intended to cover the current hybrid evaluation map:
series in broad easy `|z| < 1` regimes, direct real-axis integration for the
complete negative-real lane, and the contour backend as the general fallback.

## Usage

Report-only mode:

```bash
julia --project=perf -e 'using Pkg; Pkg.instantiate(); include("perf/speed_regression.jl")'
```

Strict regression mode:

```bash
BFL_ENFORCE_SPEED=1 julia --project=perf -e 'using Pkg; Pkg.instantiate(); include("perf/speed_regression.jl")'
```

The strict mode compares each case against a generous per-case budget in milliseconds.
Those budgets are intended to catch large accidental slowdowns rather than tiny runtime
fluctuations.

## Dense regime matrix

For broader regime coverage without turning normal tests into long-running jobs,
there is also a denser benchmarking sweep:

```bash
julia --project=perf -e 'include("perf/regime_matrix.jl")'
```

This script probes:

- a denser real complete grid for `0 < z < 1`
- incomplete cases with `b = 1e-3, 0.05, 0.5, 1.0`
- the real window `1 < z < exp(b)`
- representative complex complete and incomplete cases
- direct series-vs-contour timing comparisons for the complete case near `z = 1`

## Near-z=1 asymptotic benchmark

To compare the new complete Bose/Lerch near-`z = 1` asymptotic patch against
the previous contour fallback and a high-precision reference:

```bash
julia --project=perf -e 'include("perf/benchmark_asymptotics.jl")'
```
