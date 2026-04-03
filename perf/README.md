# Performance Regression Harness

This folder contains a small benchmark matrix for `lerch(z, s, a, b)` that probes:

- well-damped complete values
- complete values near `z = 1`
- incomplete values with small and large `b`
- the real window `1 < z < exp(b)` where the incomplete integral is valid
- complex complete and incomplete inputs

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
