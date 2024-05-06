using LumenManufaktur
using Test

@testset "direct light from muon" begin
    params = LMParameters(dispersion_model=DispersionARCA)
    pmt = LumenManufaktur.KM3NeTPMT

    # The relative tolerances are based on comparisons between
    # LumenManufaktur.jl v01.1.1 and Jpp v18.6.0-rc.1-26-g83ce39ff1
    @test isapprox(0.0304415, directlightfrommuon(params, pmt, 30, 2.33319, π), rtol=0.017)
    @test isapprox(0.00152972, directlightfrommuon(params, pmt, 23, 4.78319, π/2), rtol=0.017)
    @test isapprox(2.63515e-05, directlightfrommuon(params, pmt, 42, 1.78319, π/3), rtol=0.025)
    @test 0 ≈ directlightfrommuon(params, pmt, 5, 0.78319, 0.5708)

end
