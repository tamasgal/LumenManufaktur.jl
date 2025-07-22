using LumenManufaktur
using DelimitedFiles
using Test

@testset "direct light from muon" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = LumenManufaktur.KM3NeTPMT

    data, header = readdlm(
        "jpp/muon_direct_light_3params_Jpp_b96333ed0e136855e359ca515ea61c0138b5b43c.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = LumenManufaktur.directlightfrommuon(params, pmt, row[1], row[2], row[3])

        @test isapprox(row[4], npe, rtol=0.0001)
    end

    data, header = readdlm(
        "jpp/muon_direct_light_4params_Jpp_b96333ed0e136855e359ca515ea61c0138b5b43c.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = LumenManufaktur.directlightfrommuon(params, pmt, row[1], row[2], row[3], row[4])

        @test isapprox(row[5], npe, rtol=0.0001)
    end
end

@testset "scattered light from muon" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = LumenManufaktur.KM3NeTPMT

    data, header = readdlm(
        "jpp/muon_scattered_light_5params_Jpp_b96333ed0e136855e359ca515ea61c0138b5b43c.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfrommuon(params, pmt, row[1], row[2], row[3], row[4], row[5])
        @test isapprox(row[6], npe, atol=0.001)
    end
end
