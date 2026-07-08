# Cross-check of LumenManufaktur's exact light functions against the Jpp-generated
# PDF/CDF lookup tables (read via Lucidarium.jl).
#
# This is an opt-in integration test, NOT part of the automatic `Pkg.test()` run: it
# needs both the ~1 GB Jpp table set (`~/Data/JppPDFs`, or `KM3NET_PDF_DIR`) and
# Lucidarium.jl, neither of which is available in CI. Run it locally in an environment
# that provides both packages, e.g.
#
#   julia --project=@. test/crosscheck_jpp_tables.jl
#
# after `] dev Lucidarium`.
#
# The tables were generated for the ARCA site, so the exact calculation must use
# `DispersionARCA` and the legacy `KM3NeTPMT_Jpp` model to reproduce them. Only the
# well-behaved types are asserted here (direct light, bright point, total npe); the
# scattered secondary-source PDFs need a higher Gauss-Legendre order to converge and
# are covered qualitatively in the cross-check report, not asserted.

using LumenManufaktur
using Test

const PDF_DIR = get(ENV, "KM3NET_PDF_DIR", expanduser("~/Data/JppPDFs"))

const _HAVE_LUCIDARIUM = try
    @eval using Lucidarium
    true
catch
    false
end

const PARAMS = LMParameters(dispersion_model = DispersionARCA,
                            integration_points = PDFIntegrationPoints(5))
const PMT = KM3NeTPMT_Jpp
const FLOOR = 1.0e-7

_median(v) = (s = sort(v); n = length(s);
              n == 0 ? NaN : isodd(n) ? s[(n+1)÷2] : 0.5*(s[n÷2] + s[n÷2+1]))
_logrange(a, b, n) = exp.(range(log(a), log(b), length = n))
_pdfpath(f) = joinpath(PDF_DIR, f)

# Median ratio table/exact over the high-probability points (exact value in the top
# 20%), which is where the two implementations must agree.
function _core_ratio(pairs)
    ls = first.(pairs); ms = last.(pairs)
    isempty(ms) && return (NaN, 0)
    thr = sort(ms)[max(1, ceil(Int, 0.8 * length(ms)))]
    r = [ls[i]/ms[i] for i in eachindex(ms) if ms[i] >= thr]
    (_median(r), length(pairs))
end

if !_HAVE_LUCIDARIUM
    @warn "Lucidarium.jl not available - skipping Jpp table cross-check. Run `] dev Lucidarium`."
elseif !isfile(_pdfpath("J1p.dat.gz"))
    @warn "Jpp tables not found at $PDF_DIR - skipping cross-check. Set KM3NET_PDF_DIR."
else

@testset "direct light from muon (J1p) vs table" begin
    pdf = load_pdf(_pdfpath("J1p.dat.gz"))
    pairs = Tuple{Float64,Float64}[]
    for R in _logrange(2.0, 130.0, 7), θ in range(1.6, 3.0, length = 7),
        φ in range(0.1, 3.0, length = 4), dt in -2.0:1.0:12.0
        l = pdf(R, θ, φ, dt)
        m = directlightfrommuon(PARAMS, PMT, R, θ, φ, dt)
        (l >= FLOOR && m >= FLOOR) && push!(pairs, (l, m))
    end
    med, n = _core_ratio(pairs)
    @info "J1p direct muon" median_ratio=med npoints=n
    @test n > 500
    @test 0.85 <= med <= 1.15
    clear_pdf_cache()
end

@testset "bright point (J23p/J24p) vs table" begin
    # Direct light has a narrow time support, so it needs a finer dt window than the
    # broad scattered tail.
    for (file, fn, dtwin, lo, hi, minn) in (
        ("J23p.dat.gz", directlightfrombrightpoint, -3.0:0.5:20.0, 0.9, 1.25, 50),
        ("J24p.dat.gz", scatteredlightfrombrightpoint, -3.0:3.0:120.0, 0.9, 1.15, 200),
    )
        pdf = load_pdf(_pdfpath(file))
        pairs = Tuple{Float64,Float64}[]
        for D in _logrange(0.6, 70.0, 8), ct in range(-0.9, 0.9, length = 11), dt in dtwin
            l = pdf(D, ct, dt)
            m = fn(PARAMS, PMT, D, ct, dt)
            (l >= FLOOR && m >= FLOOR) && push!(pairs, (l, m))
        end
        med, n = _core_ratio(pairs)
        @info "bright point" file median_ratio=med npoints=n
        @test n > minn
        @test lo <= med <= hi
        clear_pdf_cache()
    end
end

if isfile(_pdfpath("I1p.dat.gz"))
    @testset "total npe: getnpe(I1p CDF) vs directlightfrommuon" begin
        cdf = load_cdf(_pdfpath("I1p.dat.gz"))
        ratios = Float64[]
        for R in _logrange(3.0, 100.0, 6), θ in range(2.0, 3.0, length = 6), φ in range(0.2, 2.8, length = 4)
            g = getnpe(cdf, R, θ, φ)
            m = directlightfrommuon(PARAMS, PMT, R, θ, φ)
            (g >= FLOOR && m >= FLOOR) && push!(ratios, g/m)
        end
        med = _median(ratios)
        @info "getnpe(I1p) vs total npe" median_ratio=med npoints=length(ratios)
        @test length(ratios) > 50
        @test 0.8 <= med <= 1.2
        clear_pdf_cache()
    end
end

end  # gating
