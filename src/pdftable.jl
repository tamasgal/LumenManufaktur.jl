# ---------------------------------------------------------------------------
# Config types
# ---------------------------------------------------------------------------

"""
    PDFTableConfig(; kwargs...)

Configuration for [`makepdf`](@ref) (`DirectMuonPDFTable`).

# Keyword Arguments
- `R_values`: grid of closest-approach distances [m]
- `θ_values`: grid of PMT zenith angles [rad]
- `ϕ_values`: grid of PMT azimuth angles [rad]
- `n_quantiles`: number of CDF quantile knots on the time axis
- `t_min`: lower bound of the time window to scan [ns]
- `t_max`: upper bound of the time window to scan [ns]
- `n_integration`: number of points for the trapezoidal CDF integration
"""
Base.@kwdef struct PDFTableConfig
    R_values::Vector{Float64} = Float64[
        0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0,
        16.0, 20.0, 25.0, 30.0, 40.0, 50.0, 65.0, 80.0, 100.0,
    ]
    θ_values::Vector{Float64} = collect(range(0.0, π; length = 13))
    ϕ_values::Vector{Float64} = collect(range(0.0, π; length = 7))
    n_quantiles::Int   = 100
    t_min::Float64     = -2.0
    t_max::Float64     = 25.0
    n_integration::Int = 5000
end

"""
    MuonPDFTableConfig(; kwargs...)

Configuration for [`makemuonpdf`](@ref) (`MuonPDFTable`).
Stores all six muon-light channels on a shared `(R, θ, ϕ)` grid.

# Keyword Arguments
- `R_values`, `θ_values`, `ϕ_values`: spatial grid
- `n_quantiles`: quantile knots per node
- `t_min`, `t_max`: time window for direct-type channels [ns]
- `t_max_scattered`: upper bound for scattered-light channels [ns]
- `n_integration`: trapezoidal integration points
"""
Base.@kwdef struct MuonPDFTableConfig
    R_values::Vector{Float64} = Float64[
        0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0,
        16.0, 20.0, 25.0, 30.0, 40.0, 50.0, 65.0, 80.0, 100.0,
    ]
    θ_values::Vector{Float64} = collect(range(0.0, π; length = 13))
    ϕ_values::Vector{Float64} = collect(range(0.0, π; length = 7))
    n_quantiles::Int     = 100
    t_min::Float64       = -2.0
    t_max::Float64       = 25.0
    t_max_scattered::Float64 = 200.0
    n_integration::Int   = 5000
end

"""
    EMShowerPDFTableConfig(; kwargs...)

Configuration for [`makeemshowerpdf`](@ref) (`EMShowerPDFTable`).
Direct and scattered EM-shower light on a `(D, cd, θ, ϕ)` grid.

# Keyword Arguments
- `D_values`: distances [m]
- `cd_values`: cosine of shower-direction/PMT angle
- `θ_values`, `ϕ_values`: PMT orientation angles [rad]
- `n_quantiles`, `t_min`, `t_max`, `t_max_scattered`, `n_integration`
"""
Base.@kwdef struct EMShowerPDFTableConfig
    D_values::Vector{Float64} = Float64[
        0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0,
        16.0, 20.0, 25.0, 30.0, 40.0, 50.0,
    ]
    cd_values::Vector{Float64} = collect(range(-1.0, 1.0; length = 9))
    θ_values::Vector{Float64}  = collect(range(0.0, π; length = 13))
    ϕ_values::Vector{Float64}  = collect(range(0.0, π; length = 7))
    n_quantiles::Int     = 100
    t_min::Float64       = -2.0
    t_max::Float64       = 25.0
    t_max_scattered::Float64 = 200.0
    n_integration::Int   = 5000
end

"""
    BrightPointPDFTableConfig(; kwargs...)

Configuration for [`makebrightpointpdf`](@ref) (`BrightPointPDFTable`).
Direct and scattered bright-point light on a `(D, ct)` grid.

# Keyword Arguments
- `D_values`: distances [m]
- `ct_values`: cosine of PMT orientation angle
- `n_quantiles`, `t_min`, `t_max`, `t_max_scattered`, `n_integration`
"""
Base.@kwdef struct BrightPointPDFTableConfig
    D_values::Vector{Float64} = Float64[
        0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0,
        16.0, 20.0, 25.0, 30.0, 40.0, 50.0,
    ]
    ct_values::Vector{Float64} = collect(range(-1.0, 1.0; length = 11))
    n_quantiles::Int     = 100
    t_min::Float64       = -2.0
    t_max::Float64       = 25.0
    t_max_scattered::Float64 = 200.0
    n_integration::Int   = 5000
end


# ---------------------------------------------------------------------------
# Table structs
# ---------------------------------------------------------------------------

"""
    DirectMuonPDFTable

Precomputed lookup table for [`directlightfrommuon`](@ref) using a CDF
reparametrisation on the time axis.

See also: [`makepdf`](@ref), [`save_pdftable`](@ref), [`load_pdftable`](@ref).
"""
struct DirectMuonPDFTable
    R_values::Vector{Float64}
    θ_values::Vector{Float64}
    ϕ_values::Vector{Float64}
    # Layout [j, iR, iθ, iϕ]: j is stride-1 for the hot-path scan
    t_knots::Array{Float64,4}
    f_vals::Array{Float64,4}
    npe::Array{Float64,3}
end

"""
    MuonPDFTable

Precomputed lookup table for all six muon-light channels on a shared
`(R, θ, ϕ)` grid:

| Field prefix | Function |
|---|---|
| `direct_muon` | [`directlightfrommuon`](@ref) |
| `scattered_muon` | [`scatteredlightfrommuon`](@ref) |
| `direct_dr` | [`directlightfromdeltarays`](@ref) |
| `scattered_dr` | [`scatteredlightfromdeltarays`](@ref) |
| `direct_em` | [`directlightfromEMshowers`](@ref) |
| `scattered_em` | [`scatteredlightfromEMshowers`](@ref) |

Delta-ray and EM-shower channels are stored per GeV (or GeV/m); multiply by
energy (or energy-loss) when evaluating.

See also: [`makemuonpdf`](@ref), [`save_muonpdftable`](@ref),
[`load_muonpdftable`](@ref).
"""
struct MuonPDFTable
    R_values::Vector{Float64}
    θ_values::Vector{Float64}
    ϕ_values::Vector{Float64}
    # Each channel: t_knots/f_vals layout [j, iR, iθ, iϕ], npe layout [iR, iθ, iϕ]
    direct_muon_t::Array{Float64,4};   direct_muon_f::Array{Float64,4};   direct_muon_npe::Array{Float64,3}
    scattered_muon_t::Array{Float64,4}; scattered_muon_f::Array{Float64,4}; scattered_muon_npe::Array{Float64,3}
    direct_dr_t::Array{Float64,4};     direct_dr_f::Array{Float64,4};     direct_dr_npe::Array{Float64,3}
    scattered_dr_t::Array{Float64,4};  scattered_dr_f::Array{Float64,4};  scattered_dr_npe::Array{Float64,3}
    direct_em_t::Array{Float64,4};     direct_em_f::Array{Float64,4};     direct_em_npe::Array{Float64,3}
    scattered_em_t::Array{Float64,4};  scattered_em_f::Array{Float64,4};  scattered_em_npe::Array{Float64,3}
end

"""
    EMShowerPDFTable

Precomputed lookup table for direct and scattered EM-shower light
([`directlightfromEMshower`](@ref) / [`scatteredlightfromEMshower`](@ref),
5-parameter per-GeV form) on a `(D, cd, θ, ϕ)` grid.

Array layout: `[j, iD, icd, iθ, iϕ]`.

See also: [`makeemshowerpdf`](@ref), [`save_emshowerpdftable`](@ref),
[`load_emshowerpdftable`](@ref).
"""
struct EMShowerPDFTable
    D_values::Vector{Float64}
    cd_values::Vector{Float64}
    θ_values::Vector{Float64}
    ϕ_values::Vector{Float64}
    # Layout [j, iD, icd, iθ, iϕ]
    direct_t::Array{Float64,5};    direct_f::Array{Float64,5};    direct_npe::Array{Float64,4}
    scattered_t::Array{Float64,5}; scattered_f::Array{Float64,5}; scattered_npe::Array{Float64,4}
end

"""
    BrightPointPDFTable

Precomputed lookup table for direct and scattered bright-point light
([`directlightfrombrightpoint`](@ref) / [`scatteredlightfrombrightpoint`](@ref))
on a `(D, ct)` grid.

Array layout: `[j, iD, ict]`.

See also: [`makebrightpointpdf`](@ref), [`save_brightpointpdftable`](@ref),
[`load_brightpointpdftable`](@ref).
"""
struct BrightPointPDFTable
    D_values::Vector{Float64}
    ct_values::Vector{Float64}
    # Layout [j, iD, ict]
    direct_t::Array{Float64,3};    direct_f::Array{Float64,3};    direct_npe::Array{Float64,2}
    scattered_t::Array{Float64,3}; scattered_f::Array{Float64,3}; scattered_npe::Array{Float64,2}
end


# ---------------------------------------------------------------------------
# show methods
# ---------------------------------------------------------------------------

function Base.show(io::IO, t::DirectMuonPDFTable)
    nR = length(t.R_values); nθ = length(t.θ_values); nϕ = length(t.ϕ_values)
    nq = size(t.t_knots, 1)
    kb = (sizeof(t.t_knots) + sizeof(t.f_vals) + sizeof(t.npe)) ÷ 1024
    print(io, "DirectMuonPDFTable($nR×$nθ×$nϕ nodes, $nq quantiles, $(kb) KiB)")
end

function Base.show(io::IO, t::MuonPDFTable)
    nR = length(t.R_values); nθ = length(t.θ_values); nϕ = length(t.ϕ_values)
    nq = size(t.direct_muon_t, 1)
    kb = (sizeof(t.direct_muon_t) + sizeof(t.direct_muon_f) + sizeof(t.direct_muon_npe) +
          sizeof(t.scattered_muon_t) + sizeof(t.scattered_muon_f) + sizeof(t.scattered_muon_npe) +
          sizeof(t.direct_dr_t) + sizeof(t.direct_dr_f) + sizeof(t.direct_dr_npe) +
          sizeof(t.scattered_dr_t) + sizeof(t.scattered_dr_f) + sizeof(t.scattered_dr_npe) +
          sizeof(t.direct_em_t) + sizeof(t.direct_em_f) + sizeof(t.direct_em_npe) +
          sizeof(t.scattered_em_t) + sizeof(t.scattered_em_f) + sizeof(t.scattered_em_npe)) ÷ 1024
    print(io, "MuonPDFTable($nR×$nθ×$nϕ nodes, $nq quantiles, 6 channels, $(kb) KiB)")
end

function Base.show(io::IO, t::EMShowerPDFTable)
    nD = length(t.D_values); ncd = length(t.cd_values)
    nθ = length(t.θ_values); nϕ = length(t.ϕ_values)
    nq = size(t.direct_t, 1)
    kb = (sizeof(t.direct_t) + sizeof(t.direct_f) + sizeof(t.direct_npe) +
          sizeof(t.scattered_t) + sizeof(t.scattered_f) + sizeof(t.scattered_npe)) ÷ 1024
    print(io, "EMShowerPDFTable($nD×$ncd×$nθ×$nϕ nodes, $nq quantiles, $(kb) KiB)")
end

function Base.show(io::IO, t::BrightPointPDFTable)
    nD = length(t.D_values); nct = length(t.ct_values)
    nq = size(t.direct_t, 1)
    kb = (sizeof(t.direct_t) + sizeof(t.direct_f) + sizeof(t.direct_npe) +
          sizeof(t.scattered_t) + sizeof(t.scattered_f) + sizeof(t.scattered_npe)) ÷ 1024
    print(io, "BrightPointPDFTable($nD×$nct nodes, $nq quantiles, $(kb) KiB)")
end


# ---------------------------------------------------------------------------
# Generic node builder: fn(Δt) -> Float64
# ---------------------------------------------------------------------------

function _build_node(fn, t_min::Float64, t_max::Float64, n_int::Int, n_q::Int)
    dt  = (t_max - t_min) / (n_int - 1)
    ts  = range(t_min, t_max; length = n_int)
    fs  = [fn(Float64(t)) for t in ts]

    cdf = zeros(n_int)
    @inbounds for i in 2:n_int
        cdf[i] = cdf[i-1] + 0.5 * (fs[i-1] + fs[i]) * dt
    end
    total_npe = cdf[end]

    t_q = Vector{Float64}(undef, n_q)
    f_q = Vector{Float64}(undef, n_q)

    if total_npe < 1e-300
        fill!(t_q, (t_min + t_max) / 2)
        fill!(f_q, 0.0)
        return t_q, f_q, 0.0
    end

    @. cdf /= total_npe

    ts_v = collect(ts)
    j = 1
    @inbounds for (k, p) in enumerate(range(0.0, 1.0; length = n_q))
        while j < n_int - 1 && cdf[j + 1] < p
            j += 1
        end
        Δc   = cdf[j + 1] - cdf[j]
        frac = Δc > 0.0 ? (p - cdf[j]) / Δc : 0.0
        t_q[k] = ts_v[j] + frac * dt
        f_q[k] = fs[j]   + frac * (fs[j + 1] - fs[j])
    end

    return t_q, f_q, total_npe
end

# Backward-compat wrapper used by makepdf (DirectMuonPDFTable)
function _build_node(params, pmt, R, θ, ϕ, t_min, t_max, n_int, n_q)
    _build_node(
        Δt -> directlightfrommuon(params, pmt, R, θ, ϕ, Δt),
        t_min, t_max, n_int, n_q,
    )
end


# ---------------------------------------------------------------------------
# Grid-location helper
# ---------------------------------------------------------------------------

@inline function _locate(vals::Vector{Float64}, v::Float64)
    n = length(vals)
    i = clamp(searchsortedlast(vals, v), 1, n - 1)
    w = clamp((v - vals[i]) / (vals[i + 1] - vals[i]), 0.0, 1.0)
    return i, w
end


# ---------------------------------------------------------------------------
# Quantile-scan evaluation for different interpolation dimensions
# ---------------------------------------------------------------------------

# 2D bilinear, array layout [nq, n1, n2]
@inline function _eval2(
    t_knots::Array{Float64,3}, f_vals::Array{Float64,3},
    Δt::Float64,
    i1::Int, i2::Int,
    w1::Float64, w2::Float64,
)
    w00 = (1-w1)*(1-w2); w10 = w1*(1-w2)
    w01 = (1-w1)*w2;     w11 = w1*w2
    @inline interp(A, j) = (
        w00*A[j,i1,  i2  ] + w10*A[j,i1+1,i2  ] +
        w01*A[j,i1,  i2+1] + w11*A[j,i1+1,i2+1]
    )
    nq = size(t_knots, 1)
    t_prev = interp(t_knots, 1); f_prev = interp(f_vals, 1)
    Δt < t_prev && return 0.0
    @inbounds for j in 2:nq
        t_curr = interp(t_knots, j); f_curr = interp(f_vals, j)
        if Δt <= t_curr
            dt   = t_curr - t_prev
            frac = dt > 0.0 ? (Δt - t_prev) / dt : 0.0
            return max(0.0, f_prev + frac * (f_curr - f_prev))
        end
        t_prev = t_curr; f_prev = f_curr
    end
    return 0.0
end

# 3D trilinear, array layout [nq, n1, n2, n3]
@inline function _eval3(
    t_knots::Array{Float64,4}, f_vals::Array{Float64,4},
    Δt::Float64,
    i1::Int, i2::Int, i3::Int,
    w1::Float64, w2::Float64, w3::Float64,
)
    w000 = (1-w1)*(1-w2)*(1-w3); w100 =    w1 *(1-w2)*(1-w3)
    w010 = (1-w1)*   w2 *(1-w3); w110 =    w1 *   w2 *(1-w3)
    w001 = (1-w1)*(1-w2)*   w3;  w101 =    w1 *(1-w2)*   w3
    w011 = (1-w1)*   w2 *   w3;  w111 =    w1 *   w2 *   w3
    @inline interp(A, j) = (
        w000*A[j,i1,  i2,  i3  ] + w100*A[j,i1+1,i2,  i3  ] +
        w010*A[j,i1,  i2+1,i3  ] + w110*A[j,i1+1,i2+1,i3  ] +
        w001*A[j,i1,  i2,  i3+1] + w101*A[j,i1+1,i2,  i3+1] +
        w011*A[j,i1,  i2+1,i3+1] + w111*A[j,i1+1,i2+1,i3+1]
    )
    nq = size(t_knots, 1)
    t_prev = interp(t_knots, 1); f_prev = interp(f_vals, 1)
    Δt < t_prev && return 0.0
    @inbounds for j in 2:nq
        t_curr = interp(t_knots, j); f_curr = interp(f_vals, j)
        if Δt <= t_curr
            dt   = t_curr - t_prev
            frac = dt > 0.0 ? (Δt - t_prev) / dt : 0.0
            return max(0.0, f_prev + frac * (f_curr - f_prev))
        end
        t_prev = t_curr; f_prev = f_curr
    end
    return 0.0
end

# 4D quadrilinear, array layout [nq, n1, n2, n3, n4]
@inline function _eval4(
    t_knots::Array{Float64,5}, f_vals::Array{Float64,5},
    Δt::Float64,
    i1::Int, i2::Int, i3::Int, i4::Int,
    w1::Float64, w2::Float64, w3::Float64, w4::Float64,
)
    a1=1-w1; a2=1-w2; a3=1-w3; a4=1-w4
    w0000=a1*a2*a3*a4; w1000=w1*a2*a3*a4; w0100=a1*w2*a3*a4; w1100=w1*w2*a3*a4
    w0010=a1*a2*w3*a4; w1010=w1*a2*w3*a4; w0110=a1*w2*w3*a4; w1110=w1*w2*w3*a4
    w0001=a1*a2*a3*w4; w1001=w1*a2*a3*w4; w0101=a1*w2*a3*w4; w1101=w1*w2*a3*w4
    w0011=a1*a2*w3*w4; w1011=w1*a2*w3*w4; w0111=a1*w2*w3*w4; w1111=w1*w2*w3*w4
    @inline interp(A, j) = (
        w0000*A[j,i1,  i2,  i3,  i4  ] + w1000*A[j,i1+1,i2,  i3,  i4  ] +
        w0100*A[j,i1,  i2+1,i3,  i4  ] + w1100*A[j,i1+1,i2+1,i3,  i4  ] +
        w0010*A[j,i1,  i2,  i3+1,i4  ] + w1010*A[j,i1+1,i2,  i3+1,i4  ] +
        w0110*A[j,i1,  i2+1,i3+1,i4  ] + w1110*A[j,i1+1,i2+1,i3+1,i4  ] +
        w0001*A[j,i1,  i2,  i3,  i4+1] + w1001*A[j,i1+1,i2,  i3,  i4+1] +
        w0101*A[j,i1,  i2+1,i3,  i4+1] + w1101*A[j,i1+1,i2+1,i3,  i4+1] +
        w0011*A[j,i1,  i2,  i3+1,i4+1] + w1011*A[j,i1+1,i2,  i3+1,i4+1] +
        w0111*A[j,i1,  i2+1,i3+1,i4+1] + w1111*A[j,i1+1,i2+1,i3+1,i4+1]
    )
    nq = size(t_knots, 1)
    t_prev = interp(t_knots, 1); f_prev = interp(f_vals, 1)
    Δt < t_prev && return 0.0
    @inbounds for j in 2:nq
        t_curr = interp(t_knots, j); f_curr = interp(f_vals, j)
        if Δt <= t_curr
            dt   = t_curr - t_prev
            frac = dt > 0.0 ? (Δt - t_prev) / dt : 0.0
            return max(0.0, f_prev + frac * (f_curr - f_prev))
        end
        t_prev = t_curr; f_prev = f_curr
    end
    return 0.0
end

# npe-only trilinear helpers (3D array)
@inline function _trilinear3(A::Array{Float64,3}, iR, iθ, iϕ, wR, wθ, wϕ)
    c00 = A[iR,  iθ,  iϕ]*(1-wR) + A[iR+1,iθ,  iϕ]*wR
    c10 = A[iR,  iθ+1,iϕ]*(1-wR) + A[iR+1,iθ+1,iϕ]*wR
    c01 = A[iR,  iθ,  iϕ+1]*(1-wR) + A[iR+1,iθ,  iϕ+1]*wR
    c11 = A[iR,  iθ+1,iϕ+1]*(1-wR) + A[iR+1,iθ+1,iϕ+1]*wR
    c0  = c00*(1-wθ) + c10*wθ
    c1  = c01*(1-wθ) + c11*wθ
    return c0*(1-wϕ) + c1*wϕ
end

@inline function _bilinear2(A::Array{Float64,2}, i1, i2, w1, w2)
    return A[i1,  i2]*(1-w1)*(1-w2) + A[i1+1,i2]*w1*(1-w2) +
           A[i1,  i2+1]*(1-w1)*w2   + A[i1+1,i2+1]*w1*w2
end


# ---------------------------------------------------------------------------
# Build helpers for multi-channel 3D tables
# ---------------------------------------------------------------------------

function _build_grid3(fns, tmax_vec, config_nq, config_t_min, config_n_int,
                      nq, nR, nθ, nϕ, Rs, θs, ϕs)
    # fns: vector of (fn(R,θ,ϕ,Δt), t_max) pairs — one per channel
    n_ch = length(fns)
    t_arr = [Array{Float64,4}(undef, nq, nR, nθ, nϕ) for _ in 1:n_ch]
    f_arr = [Array{Float64,4}(undef, nq, nR, nθ, nϕ) for _ in 1:n_ch]
    n_arr = [Array{Float64,3}(undef, nR, nθ, nϕ) for _ in 1:n_ch]

    total = nR * nθ * nϕ
    done  = Threads.Atomic{Int}(0)

    Threads.@threads for iR in 1:nR
        for iθ in 1:nθ, iϕ in 1:nϕ
            for (ch, (fn, t_max)) in enumerate(zip(fns, tmax_vec))
                tq, fq, np = _build_node(
                    Δt -> fn(Rs[iR], θs[iθ], ϕs[iϕ], Δt),
                    config_t_min, t_max, config_n_int, nq,
                )
                t_arr[ch][:, iR, iθ, iϕ] .= tq
                f_arr[ch][:, iR, iθ, iϕ] .= fq
                n_arr[ch][iR, iθ, iϕ]    = np
            end
            n = Threads.atomic_add!(done, 1) + 1
            n % 50 == 0 && @info "makepdf: $n / $total nodes done"
        end
    end
    return t_arr, f_arr, n_arr
end


# ---------------------------------------------------------------------------
# makepdf  (DirectMuonPDFTable — preserved for backward compatibility)
# ---------------------------------------------------------------------------

"""
    makepdf(params::LMParameters, pmt::PMTModel, config::PDFTableConfig = PDFTableConfig())

Build a [`DirectMuonPDFTable`](@ref) by evaluating [`directlightfrommuon`](@ref)
on the grid defined by `config`.  Uses all available Julia threads.
"""
function makepdf(
    params::LMParameters,
    pmt::PMTModel,
    config::PDFTableConfig = PDFTableConfig(),
)
    Rs = config.R_values; θs = config.θ_values; ϕs = config.ϕ_values
    nR = length(Rs); nθ = length(θs); nϕ = length(ϕs)
    nq = config.n_quantiles

    t_knots = Array{Float64,4}(undef, nq, nR, nθ, nϕ)
    f_vals  = Array{Float64,4}(undef, nq, nR, nθ, nϕ)
    npe_arr = Array{Float64,3}(undef, nR, nθ, nϕ)

    total = nR * nθ * nϕ
    done  = Threads.Atomic{Int}(0)

    Threads.@threads for iR in 1:nR
        for iθ in 1:nθ, iϕ in 1:nϕ
            tq, fq, np = _build_node(params, pmt, Rs[iR], θs[iθ], ϕs[iϕ],
                                     config.t_min, config.t_max,
                                     config.n_integration, nq)
            t_knots[:, iR, iθ, iϕ] .= tq
            f_vals[:,  iR, iθ, iϕ] .= fq
            npe_arr[iR, iθ, iϕ]    = np
            n = Threads.atomic_add!(done, 1) + 1
            n % 50 == 0 && @info "makepdf: $n / $total nodes done"
        end
    end
    @info "makepdf: complete ($total nodes)"
    DirectMuonPDFTable(Rs, θs, ϕs, t_knots, f_vals, npe_arr)
end


# ---------------------------------------------------------------------------
# makemuonpdf  (MuonPDFTable — 6 channels)
# ---------------------------------------------------------------------------

"""
    makemuonpdf(params::LMParameters, pmt::PMTModel, config::MuonPDFTableConfig = MuonPDFTableConfig())

Build a [`MuonPDFTable`](@ref) with all six muon-light channels.
Uses all available Julia threads.
"""
function makemuonpdf(
    params::LMParameters,
    pmt::PMTModel,
    config::MuonPDFTableConfig = MuonPDFTableConfig(),
)
    Rs = config.R_values; θs = config.θ_values; ϕs = config.ϕ_values
    nR = length(Rs); nθ = length(θs); nϕ = length(ϕs)
    nq = config.n_quantiles
    t0 = config.t_min; t1 = config.t_max; ts = config.t_max_scattered
    ni = config.n_integration

    fns = [
        (R,θ,ϕ,Δt) -> directlightfrommuon(params, pmt, R, θ, ϕ, Δt),
        (R,θ,ϕ,Δt) -> scatteredlightfrommuon(params, pmt, R, θ, ϕ, Δt),
        (R,θ,ϕ,Δt) -> directlightfromdeltarays(params, pmt, R, θ, ϕ, Δt),
        (R,θ,ϕ,Δt) -> scatteredlightfromdeltarays(params, pmt, R, θ, ϕ, Δt),
        (R,θ,ϕ,Δt) -> directlightfromEMshowers(params, pmt, R, θ, ϕ, Δt),
        (R,θ,ϕ,Δt) -> scatteredlightfromEMshowers(params, pmt, R, θ, ϕ, Δt),
    ]
    tmaxes = [t1, ts, t1, ts, t1, ts]

    total = nR * nθ * nϕ
    done  = Threads.Atomic{Int}(0)

    ta = [Array{Float64,4}(undef, nq, nR, nθ, nϕ) for _ in 1:6]
    fa = [Array{Float64,4}(undef, nq, nR, nθ, nϕ) for _ in 1:6]
    na = [Array{Float64,3}(undef, nR, nθ, nϕ) for _ in 1:6]

    Threads.@threads for iR in 1:nR
        for iθ in 1:nθ, iϕ in 1:nϕ
            for ch in 1:6
                tq, fq, np = _build_node(
                    Δt -> fns[ch](Rs[iR], θs[iθ], ϕs[iϕ], Δt),
                    t0, tmaxes[ch], ni, nq,
                )
                ta[ch][:, iR, iθ, iϕ] .= tq
                fa[ch][:, iR, iθ, iϕ] .= fq
                na[ch][iR, iθ, iϕ]    = np
            end
            n = Threads.atomic_add!(done, 1) + 1
            n % 50 == 0 && @info "makemuonpdf: $n / $total nodes done"
        end
    end
    @info "makemuonpdf: complete ($total nodes)"

    MuonPDFTable(
        Rs, θs, ϕs,
        ta[1], fa[1], na[1],
        ta[2], fa[2], na[2],
        ta[3], fa[3], na[3],
        ta[4], fa[4], na[4],
        ta[5], fa[5], na[5],
        ta[6], fa[6], na[6],
    )
end


# ---------------------------------------------------------------------------
# makeemshowerpdf  (EMShowerPDFTable)
# ---------------------------------------------------------------------------

"""
    makeemshowerpdf(params::LMParameters, pmt::PMTModel, config::EMShowerPDFTableConfig = EMShowerPDFTableConfig())

Build an [`EMShowerPDFTable`](@ref) for the 5-parameter (per-GeV) EM-shower PDF.
Uses all available Julia threads.
"""
function makeemshowerpdf(
    params::LMParameters,
    pmt::PMTModel,
    config::EMShowerPDFTableConfig = EMShowerPDFTableConfig(),
)
    Ds  = config.D_values; cds = config.cd_values
    θs  = config.θ_values; ϕs  = config.ϕ_values
    nD  = length(Ds); ncd = length(cds); nθ = length(θs); nϕ = length(ϕs)
    nq  = config.n_quantiles
    t0  = config.t_min; t1 = config.t_max; ts = config.t_max_scattered
    ni  = config.n_integration

    dt = Array{Float64,5}(undef, nq, nD, ncd, nθ, nϕ)
    df = Array{Float64,5}(undef, nq, nD, ncd, nθ, nϕ)
    dn = Array{Float64,4}(undef, nD, ncd, nθ, nϕ)
    st = Array{Float64,5}(undef, nq, nD, ncd, nθ, nϕ)
    sf = Array{Float64,5}(undef, nq, nD, ncd, nθ, nϕ)
    sn = Array{Float64,4}(undef, nD, ncd, nθ, nϕ)

    total = nD * ncd * nθ * nϕ
    done  = Threads.Atomic{Int}(0)

    Threads.@threads for iD in 1:nD
        for icd in 1:ncd, iθ in 1:nθ, iϕ in 1:nϕ
            tq, fq, np = _build_node(
                Δt -> directlightfromEMshower(params, pmt, Ds[iD], cds[icd], θs[iθ], ϕs[iϕ], Δt),
                t0, t1, ni, nq,
            )
            dt[:, iD, icd, iθ, iϕ] .= tq
            df[:, iD, icd, iθ, iϕ] .= fq
            dn[iD, icd, iθ, iϕ]    = np

            tq, fq, np = _build_node(
                Δt -> scatteredlightfromEMshower(params, pmt, Ds[iD], cds[icd], θs[iθ], ϕs[iϕ], Δt),
                t0, ts, ni, nq,
            )
            st[:, iD, icd, iθ, iϕ] .= tq
            sf[:, iD, icd, iθ, iϕ] .= fq
            sn[iD, icd, iθ, iϕ]    = np

            n = Threads.atomic_add!(done, 1) + 1
            n % 50 == 0 && @info "makeemshowerpdf: $n / $total nodes done"
        end
    end
    @info "makeemshowerpdf: complete ($total nodes)"
    EMShowerPDFTable(Ds, cds, θs, ϕs, dt, df, dn, st, sf, sn)
end


# ---------------------------------------------------------------------------
# makebrightpointpdf  (BrightPointPDFTable)
# ---------------------------------------------------------------------------

"""
    makebrightpointpdf(params::LMParameters, pmt::PMTModel, config::BrightPointPDFTableConfig = BrightPointPDFTableConfig())

Build a [`BrightPointPDFTable`](@ref) for direct and scattered bright-point light.
Uses all available Julia threads.
"""
function makebrightpointpdf(
    params::LMParameters,
    pmt::PMTModel,
    config::BrightPointPDFTableConfig = BrightPointPDFTableConfig(),
)
    Ds  = config.D_values; cts = config.ct_values
    nD  = length(Ds); nct = length(cts)
    nq  = config.n_quantiles
    t0  = config.t_min; t1 = config.t_max; ts = config.t_max_scattered
    ni  = config.n_integration

    dt = Array{Float64,3}(undef, nq, nD, nct)
    df = Array{Float64,3}(undef, nq, nD, nct)
    dn = Array{Float64,2}(undef, nD, nct)
    st = Array{Float64,3}(undef, nq, nD, nct)
    sf = Array{Float64,3}(undef, nq, nD, nct)
    sn = Array{Float64,2}(undef, nD, nct)

    total = nD * nct
    done  = Threads.Atomic{Int}(0)

    Threads.@threads for iD in 1:nD
        for ict in 1:nct
            tq, fq, np = _build_node(
                Δt -> directlightfrombrightpoint(params, pmt, Ds[iD], cts[ict], Δt),
                t0, t1, ni, nq,
            )
            dt[:, iD, ict] .= tq; df[:, iD, ict] .= fq; dn[iD, ict] = np

            tq, fq, np = _build_node(
                Δt -> scatteredlightfrombrightpoint(params, pmt, Ds[iD], cts[ict], Δt),
                t0, ts, ni, nq,
            )
            st[:, iD, ict] .= tq; sf[:, iD, ict] .= fq; sn[iD, ict] = np

            n = Threads.atomic_add!(done, 1) + 1
            n % 10 == 0 && @info "makebrightpointpdf: $n / $total nodes done"
        end
    end
    @info "makebrightpointpdf: complete ($total nodes)"
    BrightPointPDFTable(Ds, cts, dt, df, dn, st, sf, sn)
end


# ---------------------------------------------------------------------------
# Evaluation — DirectMuonPDFTable
# ---------------------------------------------------------------------------

"""
    (table::DirectMuonPDFTable)(R, θ, ϕ, Δt)

Evaluate the time PDF [npe/ns] at `(R, θ, ϕ, Δt)` via trilinear interpolation
and a quantile-mapped time-axis scan.  No heap allocations.
"""
function (table::DirectMuonPDFTable)(R::Float64, θ::Float64, ϕ::Float64, Δt::Float64)
    iR, wR = _locate(table.R_values, R)
    iθ, wθ = _locate(table.θ_values, θ)
    iϕ, wϕ = _locate(table.ϕ_values, ϕ)
    _eval3(table.t_knots, table.f_vals, Δt, iR, iθ, iϕ, wR, wθ, wϕ)
end

(table::DirectMuonPDFTable)(R, θ, ϕ, Δt) =
    table(Float64(R), Float64(θ), Float64(ϕ), Float64(Δt))

"""
    (table::DirectMuonPDFTable)(R, θ, ϕ)

Return total expected photo-electrons via trilinear interpolation of the `npe` array.
"""
function (table::DirectMuonPDFTable)(R::Float64, θ::Float64, ϕ::Float64)
    iR, wR = _locate(table.R_values, R)
    iθ, wθ = _locate(table.θ_values, θ)
    iϕ, wϕ = _locate(table.ϕ_values, ϕ)
    _trilinear3(table.npe, iR, iθ, iϕ, wR, wθ, wϕ)
end

(table::DirectMuonPDFTable)(R, θ, ϕ) =
    table(Float64(R), Float64(θ), Float64(ϕ))


# ---------------------------------------------------------------------------
# Evaluation — MuonPDFTable
# ---------------------------------------------------------------------------

"""
    (table::MuonPDFTable)(R, θ, ϕ, Δt)

Return a `NamedTuple` with the six PDF values [npe/ns] at `(R, θ, ϕ, Δt)`.
Delta-ray and EM-shower channels are per GeV/m and per GeV respectively.
"""
function (table::MuonPDFTable)(R::Float64, θ::Float64, ϕ::Float64, Δt::Float64)
    iR, wR = _locate(table.R_values, R)
    iθ, wθ = _locate(table.θ_values, θ)
    iϕ, wϕ = _locate(table.ϕ_values, ϕ)
    (
        direct_muon      = _eval3(table.direct_muon_t,    table.direct_muon_f,    Δt, iR,iθ,iϕ, wR,wθ,wϕ),
        scattered_muon   = _eval3(table.scattered_muon_t, table.scattered_muon_f, Δt, iR,iθ,iϕ, wR,wθ,wϕ),
        direct_dr        = _eval3(table.direct_dr_t,      table.direct_dr_f,      Δt, iR,iθ,iϕ, wR,wθ,wϕ),
        scattered_dr     = _eval3(table.scattered_dr_t,   table.scattered_dr_f,   Δt, iR,iθ,iϕ, wR,wθ,wϕ),
        direct_em        = _eval3(table.direct_em_t,      table.direct_em_f,      Δt, iR,iθ,iϕ, wR,wθ,wϕ),
        scattered_em     = _eval3(table.scattered_em_t,   table.scattered_em_f,   Δt, iR,iθ,iϕ, wR,wθ,wϕ),
    )
end

(table::MuonPDFTable)(R, θ, ϕ, Δt) =
    table(Float64(R), Float64(θ), Float64(ϕ), Float64(Δt))

"""
    (table::MuonPDFTable)(R, θ, ϕ, E, Δt)

Return total muon-light PDF [npe/ns] for muon energy `E` [GeV], combining all
six channels with the same weighting as [`lightfrommuon`](@ref).
"""
function (table::MuonPDFTable)(R::Float64, θ::Float64, ϕ::Float64, E::Float64, Δt::Float64)
    c = table(R, θ, ϕ, Δt)
    dE = deltarayenergyloss(E)
    return (c.direct_muon + c.scattered_muon +
            E  * (c.direct_em    + c.scattered_em) +
            dE * (c.direct_dr    + c.scattered_dr))
end

(table::MuonPDFTable)(R, θ, ϕ, E, Δt) =
    table(Float64(R), Float64(θ), Float64(ϕ), Float64(E), Float64(Δt))

"""
    (table::MuonPDFTable)(R, θ, ϕ)

Return a `NamedTuple` with the six integrated npe values at `(R, θ, ϕ)`.
"""
function (table::MuonPDFTable)(R::Float64, θ::Float64, ϕ::Float64)
    iR, wR = _locate(table.R_values, R)
    iθ, wθ = _locate(table.θ_values, θ)
    iϕ, wϕ = _locate(table.ϕ_values, ϕ)
    (
        direct_muon    = _trilinear3(table.direct_muon_npe,    iR,iθ,iϕ, wR,wθ,wϕ),
        scattered_muon = _trilinear3(table.scattered_muon_npe, iR,iθ,iϕ, wR,wθ,wϕ),
        direct_dr      = _trilinear3(table.direct_dr_npe,      iR,iθ,iϕ, wR,wθ,wϕ),
        scattered_dr   = _trilinear3(table.scattered_dr_npe,   iR,iθ,iϕ, wR,wθ,wϕ),
        direct_em      = _trilinear3(table.direct_em_npe,      iR,iθ,iϕ, wR,wθ,wϕ),
        scattered_em   = _trilinear3(table.scattered_em_npe,   iR,iθ,iϕ, wR,wθ,wϕ),
    )
end

(table::MuonPDFTable)(R, θ, ϕ) =
    table(Float64(R), Float64(θ), Float64(ϕ))


# ---------------------------------------------------------------------------
# Evaluation — EMShowerPDFTable
# ---------------------------------------------------------------------------

"""
    (table::EMShowerPDFTable)(D, cd, θ, ϕ, Δt)

Return a `NamedTuple` `(direct, scattered)` with the per-GeV PDF values
[npe/ns/GeV] at `(D, cd, θ, ϕ, Δt)`.
"""
function (table::EMShowerPDFTable)(D::Float64, cd::Float64, θ::Float64, ϕ::Float64, Δt::Float64)
    iD,  wD  = _locate(table.D_values,  D)
    icd, wcd = _locate(table.cd_values, cd)
    iθ,  wθ  = _locate(table.θ_values,  θ)
    iϕ,  wϕ  = _locate(table.ϕ_values,  ϕ)
    (
        direct    = _eval4(table.direct_t,    table.direct_f,    Δt, iD,icd,iθ,iϕ, wD,wcd,wθ,wϕ),
        scattered = _eval4(table.scattered_t, table.scattered_f, Δt, iD,icd,iθ,iϕ, wD,wcd,wθ,wϕ),
    )
end

(table::EMShowerPDFTable)(D, cd, θ, ϕ, Δt) =
    table(Float64(D), Float64(cd), Float64(θ), Float64(ϕ), Float64(Δt))

"""
    (table::EMShowerPDFTable)(D, cd, θ, ϕ, E, Δt)

Return the total EM-shower PDF [npe/ns] scaled by shower energy `E` [GeV].
"""
function (table::EMShowerPDFTable)(D::Float64, cd::Float64, θ::Float64, ϕ::Float64, E::Float64, Δt::Float64)
    c = table(D, cd, θ, ϕ, Δt)
    return E * (c.direct + c.scattered)
end

(table::EMShowerPDFTable)(D, cd, θ, ϕ, E, Δt) =
    table(Float64(D), Float64(cd), Float64(θ), Float64(ϕ), Float64(E), Float64(Δt))

"""
    (table::EMShowerPDFTable)(D, cd, θ, ϕ)

Return `(direct, scattered)` integrated npe [per GeV] at `(D, cd, θ, ϕ)`.
"""
function (table::EMShowerPDFTable)(D::Float64, cd::Float64, θ::Float64, ϕ::Float64)
    iD,  wD  = _locate(table.D_values,  D)
    icd, wcd = _locate(table.cd_values, cd)
    iθ,  wθ  = _locate(table.θ_values,  θ)
    iϕ,  wϕ  = _locate(table.ϕ_values,  ϕ)
    a1=1-wD; a2=1-wcd; a3=1-wθ; a4=1-wϕ
    @inline quad4(A) = (
        A[iD,  icd,  iθ,  iϕ  ]*a1*a2*a3*a4 + A[iD+1,icd,  iθ,  iϕ  ]*wD*a2*a3*a4 +
        A[iD,  icd+1,iθ,  iϕ  ]*a1*wcd*a3*a4 + A[iD+1,icd+1,iθ,  iϕ  ]*wD*wcd*a3*a4 +
        A[iD,  icd,  iθ+1,iϕ  ]*a1*a2*wθ*a4 + A[iD+1,icd,  iθ+1,iϕ  ]*wD*a2*wθ*a4 +
        A[iD,  icd+1,iθ+1,iϕ  ]*a1*wcd*wθ*a4 + A[iD+1,icd+1,iθ+1,iϕ  ]*wD*wcd*wθ*a4 +
        A[iD,  icd,  iθ,  iϕ+1]*a1*a2*a3*wϕ + A[iD+1,icd,  iθ,  iϕ+1]*wD*a2*a3*wϕ +
        A[iD,  icd+1,iθ,  iϕ+1]*a1*wcd*a3*wϕ + A[iD+1,icd+1,iθ,  iϕ+1]*wD*wcd*a3*wϕ +
        A[iD,  icd,  iθ+1,iϕ+1]*a1*a2*wθ*wϕ + A[iD+1,icd,  iθ+1,iϕ+1]*wD*a2*wθ*wϕ +
        A[iD,  icd+1,iθ+1,iϕ+1]*a1*wcd*wθ*wϕ + A[iD+1,icd+1,iθ+1,iϕ+1]*wD*wcd*wθ*wϕ
    )
    (direct = quad4(table.direct_npe), scattered = quad4(table.scattered_npe))
end

(table::EMShowerPDFTable)(D, cd, θ, ϕ) =
    table(Float64(D), Float64(cd), Float64(θ), Float64(ϕ))


# ---------------------------------------------------------------------------
# Evaluation — BrightPointPDFTable
# ---------------------------------------------------------------------------

"""
    (table::BrightPointPDFTable)(D, ct, Δt)

Return a `NamedTuple` `(direct, scattered)` with the PDF values [npe/ns]
at `(D, ct, Δt)`.
"""
function (table::BrightPointPDFTable)(D::Float64, ct::Float64, Δt::Float64)
    iD,  wD  = _locate(table.D_values,  D)
    ict, wct = _locate(table.ct_values, ct)
    (
        direct    = _eval2(table.direct_t,    table.direct_f,    Δt, iD, ict, wD, wct),
        scattered = _eval2(table.scattered_t, table.scattered_f, Δt, iD, ict, wD, wct),
    )
end

(table::BrightPointPDFTable)(D, ct, Δt) =
    table(Float64(D), Float64(ct), Float64(Δt))

"""
    (table::BrightPointPDFTable)(D, ct)

Return `(direct, scattered)` integrated npe at `(D, ct)`.
"""
function (table::BrightPointPDFTable)(D::Float64, ct::Float64)
    iD,  wD  = _locate(table.D_values,  D)
    ict, wct = _locate(table.ct_values, ct)
    (
        direct    = _bilinear2(table.direct_npe,    iD, ict, wD, wct),
        scattered = _bilinear2(table.scattered_npe, iD, ict, wD, wct),
    )
end

(table::BrightPointPDFTable)(D, ct) = table(Float64(D), Float64(ct))


# ---------------------------------------------------------------------------
# HDF5 persistence — DirectMuonPDFTable
# ---------------------------------------------------------------------------

"""
    save_pdftable(filename, table)

Write a [`DirectMuonPDFTable`](@ref) to an HDF5 file.
"""
function save_pdftable(filename::AbstractString, table::DirectMuonPDFTable)
    HDF5.h5open(filename, "w") do fid
        fid["R_values"]     = table.R_values
        fid["theta_values"] = table.θ_values
        fid["phi_values"]   = table.ϕ_values
        fid["t_knots"]      = table.t_knots
        fid["f_vals"]       = table.f_vals
        fid["npe"]          = table.npe
        HDF5.attributes(fid)["type"]    = "DirectMuonPDFTable"
        HDF5.attributes(fid)["version"] = 1
    end
    return filename
end

"""
    load_pdftable(filename)

Read a [`DirectMuonPDFTable`](@ref) from an HDF5 file written by [`save_pdftable`](@ref).
"""
function load_pdftable(filename::AbstractString)
    HDF5.h5open(filename, "r") do fid
        DirectMuonPDFTable(
            read(fid["R_values"]), read(fid["theta_values"]), read(fid["phi_values"]),
            read(fid["t_knots"]), read(fid["f_vals"]), read(fid["npe"]),
        )
    end
end


# ---------------------------------------------------------------------------
# HDF5 persistence — MuonPDFTable
# ---------------------------------------------------------------------------

"""
    save_muonpdftable(filename, table)

Write a [`MuonPDFTable`](@ref) to an HDF5 file.
"""
function save_muonpdftable(filename::AbstractString, table::MuonPDFTable)
    HDF5.h5open(filename, "w") do fid
        fid["R_values"]     = table.R_values
        fid["theta_values"] = table.θ_values
        fid["phi_values"]   = table.ϕ_values
        for (prefix, t, f, n) in (
            ("direct_muon",    table.direct_muon_t,    table.direct_muon_f,    table.direct_muon_npe),
            ("scattered_muon", table.scattered_muon_t, table.scattered_muon_f, table.scattered_muon_npe),
            ("direct_dr",      table.direct_dr_t,      table.direct_dr_f,      table.direct_dr_npe),
            ("scattered_dr",   table.scattered_dr_t,   table.scattered_dr_f,   table.scattered_dr_npe),
            ("direct_em",      table.direct_em_t,      table.direct_em_f,      table.direct_em_npe),
            ("scattered_em",   table.scattered_em_t,   table.scattered_em_f,   table.scattered_em_npe),
        )
            fid["$(prefix)_t"] = t
            fid["$(prefix)_f"] = f
            fid["$(prefix)_npe"] = n
        end
        HDF5.attributes(fid)["type"]    = "MuonPDFTable"
        HDF5.attributes(fid)["version"] = 1
    end
    return filename
end

"""
    load_muonpdftable(filename)

Read a [`MuonPDFTable`](@ref) from an HDF5 file written by [`save_muonpdftable`](@ref).
"""
function load_muonpdftable(filename::AbstractString)
    HDF5.h5open(filename, "r") do fid
        MuonPDFTable(
            read(fid["R_values"]), read(fid["theta_values"]), read(fid["phi_values"]),
            read(fid["direct_muon_t"]),    read(fid["direct_muon_f"]),    read(fid["direct_muon_npe"]),
            read(fid["scattered_muon_t"]), read(fid["scattered_muon_f"]), read(fid["scattered_muon_npe"]),
            read(fid["direct_dr_t"]),      read(fid["direct_dr_f"]),      read(fid["direct_dr_npe"]),
            read(fid["scattered_dr_t"]),   read(fid["scattered_dr_f"]),   read(fid["scattered_dr_npe"]),
            read(fid["direct_em_t"]),      read(fid["direct_em_f"]),      read(fid["direct_em_npe"]),
            read(fid["scattered_em_t"]),   read(fid["scattered_em_f"]),   read(fid["scattered_em_npe"]),
        )
    end
end


# ---------------------------------------------------------------------------
# HDF5 persistence — EMShowerPDFTable
# ---------------------------------------------------------------------------

"""
    save_emshowerpdftable(filename, table)

Write an [`EMShowerPDFTable`](@ref) to an HDF5 file.
"""
function save_emshowerpdftable(filename::AbstractString, table::EMShowerPDFTable)
    HDF5.h5open(filename, "w") do fid
        fid["D_values"]     = table.D_values
        fid["cd_values"]    = table.cd_values
        fid["theta_values"] = table.θ_values
        fid["phi_values"]   = table.ϕ_values
        fid["direct_t"]   = table.direct_t;   fid["direct_f"]   = table.direct_f;   fid["direct_npe"]   = table.direct_npe
        fid["scattered_t"] = table.scattered_t; fid["scattered_f"] = table.scattered_f; fid["scattered_npe"] = table.scattered_npe
        HDF5.attributes(fid)["type"]    = "EMShowerPDFTable"
        HDF5.attributes(fid)["version"] = 1
    end
    return filename
end

"""
    load_emshowerpdftable(filename)

Read an [`EMShowerPDFTable`](@ref) from an HDF5 file.
"""
function load_emshowerpdftable(filename::AbstractString)
    HDF5.h5open(filename, "r") do fid
        EMShowerPDFTable(
            read(fid["D_values"]), read(fid["cd_values"]),
            read(fid["theta_values"]), read(fid["phi_values"]),
            read(fid["direct_t"]),   read(fid["direct_f"]),   read(fid["direct_npe"]),
            read(fid["scattered_t"]), read(fid["scattered_f"]), read(fid["scattered_npe"]),
        )
    end
end


# ---------------------------------------------------------------------------
# HDF5 persistence — BrightPointPDFTable
# ---------------------------------------------------------------------------

"""
    save_brightpointpdftable(filename, table)

Write a [`BrightPointPDFTable`](@ref) to an HDF5 file.
"""
function save_brightpointpdftable(filename::AbstractString, table::BrightPointPDFTable)
    HDF5.h5open(filename, "w") do fid
        fid["D_values"]  = table.D_values
        fid["ct_values"] = table.ct_values
        fid["direct_t"]   = table.direct_t;   fid["direct_f"]   = table.direct_f;   fid["direct_npe"]   = table.direct_npe
        fid["scattered_t"] = table.scattered_t; fid["scattered_f"] = table.scattered_f; fid["scattered_npe"] = table.scattered_npe
        HDF5.attributes(fid)["type"]    = "BrightPointPDFTable"
        HDF5.attributes(fid)["version"] = 1
    end
    return filename
end

"""
    load_brightpointpdftable(filename)

Read a [`BrightPointPDFTable`](@ref) from an HDF5 file.
"""
function load_brightpointpdftable(filename::AbstractString)
    HDF5.h5open(filename, "r") do fid
        BrightPointPDFTable(
            read(fid["D_values"]), read(fid["ct_values"]),
            read(fid["direct_t"]),   read(fid["direct_f"]),   read(fid["direct_npe"]),
            read(fid["scattered_t"]), read(fid["scattered_f"]), read(fid["scattered_npe"]),
        )
    end
end
