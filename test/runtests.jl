using LumenManufaktur
using DelimitedFiles
using BasicInterpolators
using Test


# KM3NeTPMT_Jpp is exported by LumenManufaktur (see src/pmt.jl).



@testset "direct light from muon" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/muon_direct_light_3params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = LumenManufaktur.directlightfrommuon(params, pmt, row[1], row[2], row[3])

        @test isapprox(row[4], npe, rtol=0.0001)
    end

    data, header = readdlm(
        "jpp/muon_direct_light_4params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = LumenManufaktur.directlightfrommuon(params, pmt, row[1], row[2], row[3], row[4])

        @test isapprox(row[5], npe, rtol=0.0001)
    end
end

@testset "scattered light from muon (5-param, per-segment)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/muon_scattered_light_5params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfrommuon(params, pmt, row[1], row[2], row[3], row[4], row[5])
        @test isapprox(row[6], npe, atol=0.001)
    end
end

@testset "scattered light from muon (4-param, track-integrated)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/muon_scattered_light_4params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfrommuon(params, pmt, row[1], row[2], row[3], row[4])
        @test isapprox(row[5], npe, atol=0.001)
    end
end

@testset "direct light from delta rays" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/delta_rays_direct_light_4params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = directlightfromdeltarays(params, pmt, row[1], row[2], row[3], row[4])
        @test isapprox(row[5], npe, rtol=0.0001)
    end
end

@testset "scattered light from delta rays" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/delta_rays_scattered_light_4params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfromdeltarays(params, pmt, row[1], row[2], row[3], row[4])
        @test isapprox(row[5], npe, atol=0.001)
    end
end

@testset "direct light from EM-shower (5-param)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/emshower_direct_light_5params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = directlightfromEMshower(params, pmt, row[1], row[2], row[3], row[4], row[5])
        # ~0.12% systematic offset from _lgamma approximation vs C++ tgamma; rtol=0.002 accommodates this
        @test isapprox(row[6], npe, rtol=0.002)
    end
end

@testset "scattered light from EM-shower (5-param)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/emshower_scattered_light_5params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfromEMshower(params, pmt, row[1], row[2], row[3], row[4], row[5])
        @test isapprox(row[6], npe, atol=0.001)
    end
end

@testset "direct light from EM-shower (6-param, with energy)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/emshower_direct_light_6params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = directlightfromEMshower(params, pmt, row[1], row[2], row[3], row[4], row[5], row[6])
        # ~0.33% max offset: the 6-param function integrates 5-param calls, each with ~0.12%
        # systematic offset (C++ uses tabulated geant vs Julia's exact evaluation), accumulating
        # to up to 0.33%; rtol=0.004 accommodates this
        @test isapprox(row[7], npe, rtol=0.004)
    end
end

@testset "scattered light from EM-shower (6-param, with energy)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/emshower_scattered_light_6params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfromEMshower(params, pmt, row[1], row[2], row[3], row[4], row[5], row[6])
        @test isapprox(row[7], npe, atol=0.001)
    end
end

@testset "direct light from EM-showers (track-integrated)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/emshowers_direct_light_4params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = directlightfromEMshowers(params, pmt, row[1], row[2], row[3], row[4])
        @test isapprox(row[5], npe, rtol=0.002)
    end
end

@testset "scattered light from EM-showers (track-integrated)" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/emshowers_scattered_light_4params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfromEMshowers(params, pmt, row[1], row[2], row[3], row[4])
        @test isapprox(row[5], npe, atol=0.001)
    end
end

@testset "direct light from bright point" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/brightpoint_direct_light_3params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = directlightfrombrightpoint(params, pmt, row[1], row[2], row[3])
        @test isapprox(row[4], npe, rtol=0.0001)
    end
end

@testset "scattered light from bright point" begin
    params = LMParameters(
        dispersion_model = DispersionORCA,
        integration_points = PDFIntegrationPoints(5),
    )
    pmt = KM3NeTPMT_Jpp

    data, header = readdlm(
        "jpp/brightpoint_scattered_light_3params_Jpp_ad6bf9c8a.csv",
        ',',
        Float64,
        '\n'; header=true
    )
    for row in eachrow(data)
        npe = scatteredlightfrombrightpoint(params, pmt, row[1], row[2], row[3])
        @test isapprox(row[4], npe, atol=0.001)
    end
end
