using Test
using QuadGK
using SpecialFunctions

# This file validates the manuscript-style near-degeneracy Bose claims using
# the public API and an independent high-precision quadrature reference.

function big_bose_quad_ref(z, s; rtol = big"1e-28", prec = 256)
    @assert 0 < z < 1
    @assert s > 0
    setprecision(BigFloat, prec) do
        zb = BigFloat(z)
        sb = BigFloat(s)
        gb = gamma(sb)
        integrand(t) = zb * exp(-t) * t^(sb - 1) / ((1 - zb * exp(-t)) * gb)
        quadgk(integrand, big"0.0", big"1.0", big"4.0", big"16.0", big"64.0", Inf; rtol = rtol)[1]
    end
end

relerr(x, y) = abs(x - y) / abs(y)

maybe_public_bose(z, s) = bose(z, s; rtol = 1e-12)

@testset "bose near z=1 public accuracy" begin
    noninteger_orders = (0.5, 1.5, 2.5)
    integer_orders = (1.0, 2.0, 3.0)
    near_one_z = (0.99, 0.999, 0.9999, 0.99999)

    @testset "noninteger accuracy near z=1" begin
        for s in noninteger_orders, z in near_one_z
            y_pkg = maybe_public_bose(z, s)
            y_ref = Float64(big_bose_quad_ref(z, s))

            @test isfinite(y_pkg)
            @test y_pkg > 0
            @test relerr(y_pkg, y_ref) <= 1e-10
        end
    end

    @testset "integer accuracy near z=1" begin
        for z in near_one_z
            y_pkg = maybe_public_bose(z, 1.0)
            y_ref = -log1p(-z)

            @test isfinite(y_pkg)
            @test y_pkg > 0
            @test relerr(y_pkg, y_ref) <= 1e-13
        end

        for s in (2.0, 3.0), z in near_one_z
            y_pkg = maybe_public_bose(z, s)
            y_ref = Float64(big_bose_quad_ref(z, s))

            @test isfinite(y_pkg)
            @test y_pkg > 0
            @test relerr(y_pkg, y_ref) <= 1e-10
        end
    end

    @testset "switching-region continuity" begin
        for s in (0.5, 1.5)
            for z in (0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
                y_pkg = maybe_public_bose(z, s)
                y_ref = Float64(big_bose_quad_ref(z, s))

                @test isfinite(y_pkg)
                @test y_pkg > 0
                @test relerr(y_pkg, y_ref) <= 5e-10
            end
        end
    end

    @testset "monotonicity in z" begin
        zgrid = (0.9, 0.99, 0.999, 0.9999)
        for s in (noninteger_orders..., integer_orders...)
            vals = [maybe_public_bose(z, s) for z in zgrid]
            @test all(isfinite, vals)
            @test all(>(0), vals)
            @test all(diff(vals) .> 0)
        end
    end

    @testset "optional internal asymptotic consistency" begin
        if isdefined(BoseFermiLerch, :_bose_asymp_z1_real)
            for s in noninteger_orders, z in (0.99, 0.9999, 0.99999)
                y_pkg = maybe_public_bose(z, s)
                y_asymp = BoseFermiLerch._bose_asymp_z1_real(z, s; rtol = 1e-12)
                @test isapprox(y_pkg, y_asymp; rtol = 1e-12, atol = 0.0)
            end
        else
            @test true
        end
    end
end
