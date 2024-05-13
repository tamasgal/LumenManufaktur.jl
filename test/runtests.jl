using LumenManufaktur
using Test

@testset "direct light from muon" begin
    params = LMParameters(dispersion_model = DispersionORCA)
    pmt = LumenManufaktur.KM3NeTPMT

    # The relative tolerances are based on comparisons between
    # LumenManufaktur.jl v0.5.1 and Jpp v18.6.0-rc.1-26-g83ce39ff1
    @test isapprox(
        0.0304415,
        directlightfrommuon(params, pmt, 30, 2.33319, π),
        rtol = 0.000001,
    )
    @test isapprox(
        0.00152972,
        directlightfrommuon(params, pmt, 23, 4.78319, π / 2),
        rtol = 0.0001,
    )
    @test isapprox(
        2.63515e-05,
        directlightfrommuon(params, pmt, 42, 1.78319, π / 3),
        rtol = 0.0001,
    )
    @test isapprox(
        13.0251,
        directlightfrommuon(params, pmt, 1.15, 2.85599, π / 2, 0.01),
        rtol = 0.000001,
    )
    @test 0 ≈ directlightfrommuon(params, pmt, 5, 0.78319, 0.5708)


    @test isapprox(
        0.00414637,
        scatteredlightfrommuon(params, pmt, 2.0, 0.5, π / 2, π, 1.23),
        rtol = 0.000001,
    )
#
    @test isapprox(
        0.00432646,
        directlightfromdeltarays(params, pmt, 1.5, π/2, π/3, 1.23),
        rtol = 0.000001
    )

end
