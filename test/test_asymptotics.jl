using Test
using QuadGK

fd_ref_asymp(s, z; rtol = 1e-8) = z * quadgk(t -> exp(-t) * t^(s - 1) / ((1 + z * exp(-t)) * gamma(s)), 0.0, Inf; rtol = rtol)[1]

function big_bose_ref(z, s; rtol = big"1e-30")
    setprecision(256) do
        zb = BigFloat(z)
        sb = BigFloat(s)
        gb = gamma(sb)
        integrand(t) = exp(-t) * t^(sb - 1) * zb / ((1 - zb * exp(-t)) * gb)
        quadgk(integrand, big"0.0", big"1.0", big"4.0", big"16.0", big"64.0", Inf; rtol = rtol)[1]
    end
end

function big_lerch_ref(z, s; rtol = big"1e-30")
    setprecision(256) do
        zb = BigFloat(z)
        sb = BigFloat(s)
        gb = gamma(sb)
        integrand(t) = exp(-t) * t^(sb - 1) / ((1 - zb * exp(-t)) * gb)
        quadgk(integrand, big"0.0", big"1.0", big"4.0", big"16.0", big"64.0", Inf; rtol = rtol)[1]
    end
end

@testset "near-z=1 asymptotic patch" begin
    noninteger_orders = (0.5, 1.5, 2.5, 3.7)
    near_one_z = (0.98, 0.999, 0.9999, 0.99999)

    @testset "noninteger Bose/Lerch accuracy" begin
        for s in noninteger_orders, z in near_one_z
            bose_ref = Float64(big_bose_ref(z, s))
            lerch_ref = Float64(big_lerch_ref(z, s))

            @test isapprox(BoseFermiLerch._bose_asymp_z1_real(z, s; rtol = 1e-13), bose_ref; rtol = 5e-12)
            @test isapprox(BoseFermiLerch._lerch_complete_asymp_z1_real(z, s; rtol = 1e-13), lerch_ref; rtol = 5e-12)
            @test isapprox(bose(z, s; rtol = 1e-12), bose_ref; rtol = 5e-12)
            @test isapprox(lerch(z, s, 1.0, 0.0; rtol = 1e-12), lerch_ref; rtol = 5e-12)
        end
    end

    @testset "integer-order handling" begin
        for z in near_one_z
            @test isapprox(BoseFermiLerch._bose_asymp_z1_real(z, 1.0; rtol = 1e-13), -log1p(-z); rtol = 1e-14)
            @test isapprox(BoseFermiLerch._lerch_complete_asymp_z1_real(z, 1.0; rtol = 1e-13), -log1p(-z) / z; rtol = 1e-14)
        end

        for s in (2.0, 3.0), z in near_one_z
            bose_ref = Float64(big_bose_ref(z, s))
            lerch_ref = Float64(big_lerch_ref(z, s))

            @test isapprox(BoseFermiLerch._bose_asymp_z1_real(z, s; rtol = 1e-13), bose_ref; rtol = 5e-12)
            @test isapprox(BoseFermiLerch._lerch_complete_asymp_z1_real(z, s; rtol = 1e-13), lerch_ref; rtol = 5e-12)
        end
    end

    @testset "switch continuity near asymptotic threshold" begin
        for s in (1.5, 2.5)
            z = exp(-0.05)
            contour_bose = z * BoseFermiLerch.lerch_complete_contour(z, s, 1.0; rtol = 1e-12)
            asymp_bose = BoseFermiLerch._bose_asymp_z1_real(z, s; rtol = 1e-12)
            contour_lerch = BoseFermiLerch.lerch_complete_contour(z, s, 1.0; rtol = 1e-12)
            asymp_lerch = BoseFermiLerch._lerch_complete_asymp_z1_real(z, s; rtol = 1e-12)

            @test isapprox(asymp_bose, contour_bose; rtol = 5e-10)
            @test isapprox(asymp_lerch, contour_lerch; rtol = 5e-10)
        end
    end

    @testset "common Bose orders near unity" begin
        z = 0.999999
        for s in (0.5, 1.5, 2.5)
            ref = Float64(big_bose_ref(z, s))
            @test isapprox(bose(z, s; rtol = 1e-12), ref; rtol = 5e-11)
            @test bose(z, s; rtol = 1e-12) == BoseFermiLerch._bose_asymp_z1_real(z, s; rtol = 1e-12)
        end
    end

    @testset "no regression outside patch region" begin
        @test isapprox(bose(0.8, 1.5; rtol = 1e-12), Float64(big_bose_ref(0.8, 1.5)); rtol = 1e-11)
        @test isapprox(lerch(0.8, 1.5, 1.0, 0.0; rtol = 1e-12), Float64(big_lerch_ref(0.8, 1.5)); rtol = 1e-11)
        @test isapprox(fermi(exp(30), 3 / 2; rtol = 1e-7), fd_ref_asymp(3 / 2, exp(30); rtol = 1e-7); rtol = 1e-7)
    end
end
