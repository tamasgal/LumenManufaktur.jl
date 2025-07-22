using LumenManufaktur
using DelimitedFiles
using Test

@testset "direct light from muon" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = LumenManufaktur.KM3NeTPMT

    # The relative tolerances are based on comparisons between
    # LumenManufaktur.jl v0.5.1 and Jpp v18.6.0-rc.1-26-g83ce39ff1
    @test isapprox(
        0.0339756,
        directlightfrommuon(params, pmt, 30, 2.33319, π),
        rtol = 0.00001,
    )
    @test isapprox(
        0.00170725,
        directlightfrommuon(params, pmt, 23, 4.78319, π / 2),
        rtol = 0.0001,
    )
    @test isapprox(
        2.94122e-05,
        directlightfrommuon(params, pmt, 42, 1.78319, π / 3),
        rtol = 0.0001,
    )
    @test isapprox(
        14.5372,
        directlightfrommuon(params, pmt, 1.15, 2.85599, π / 2, 0.01),
        rtol = 0.00001,
    )
    @test 0 ≈ directlightfrommuon(params, pmt, 5, 0.78319, 0.5708)


    @test isapprox(
        0.00462774,
        scatteredlightfrommuon(params, pmt, 2.0, 0.5, π / 2, π, 1.23),
        rtol = 0.000001,
    )

    @test isapprox(
        0.00482852,
        directlightfromdeltarays(params, pmt, 1.5, π / 2, π / 3, 1.23),
        rtol = 0.000001,
    )
    @test isapprox(
        0.0031653,
        scatteredlightfromdeltarays(params, pmt, 1.5, π / 2, π / 3, 1.23),
        rtol = 0.00001,
    )
    # @test isapprox(
    #     0.000394462,
    #     scatteredlightfromdeltarays(params, pmt, 2.5, π / 3, π / 4, 2.46),
    #     rtol = 0.000001,
    # )
end

@testset "brightpoint" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = LumenManufaktur.KM3NeTPMT

    data, header = readdlm(
        "brightpoint_Jpp_b96333ed0e136855e359ca515ea61c0138b5b43c.csv",
        ',',
        Float64,
        '\n'; header=true
    )


    for row in eachrow(data)
        dl = LumenManufaktur.directlightfrombrightpoint(params, pmt, row[1], row[2], row[3])
        sl = LumenManufaktur.scatteredlightfrombrightpoint(params, pmt, row[1], row[2], row[3])
        @test isapprox(row[4], dl, rtol=0.0001)
        @test isapprox(row[5], sl, rtol=0.001)
    end
end
