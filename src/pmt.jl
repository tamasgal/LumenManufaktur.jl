"""
The PMT model.

# Arguments
- `photocathode_area`: the total are of the photocathode [m^2]
- `quantum_efficiency`: a callable which returns the quantum efficiency
  for a given wavelength [ns]
- `angular_acceptance`: a callable which returns the angular acceptance
  (for values between -1 and 1) for a given wavelength [ns]
"""
struct PMTModel{T1,T2}
    photocathode_area::Float64
    quantum_efficiency::T1
    angular_acceptance::T2
end

const KM3NeTPMT = PMTModel(
    45.4e-4,
    # The quantum efficiency includes absorption in glass and gel.
    # Yields slightly different values compared to the Jpp getQE() function,
    # very probably due to the interpolation implementation.
    LinearInterpolator(
        [
            0,
            270,
            275,
            280,
            285,
            290,
            295,
            300,
            305,
            310,
            315,
            320,
            325,
            330,
            335,
            340,
            345,
            350,
            355,
            360,
            365,
            370,
            375,
            380,
            385,
            390,
            395,
            400,
            405,
            410,
            415,
            420,
            425,
            430,
            435,
            440,
            445,
            450,
            455,
            460,
            465,
            470,
            475,
            480,
            485,
            490,
            495,
            500,
            505,
            510,
            515,
            520,
            525,
            530,
            535,
            540,
            545,
            550,
            555,
            560,
            565,
            570,
            575,
            580,
            585,
            590,
            595,
            600,
            605,
            610,
            615,
            620,
            625,
            630,
            635,
            640,
            645,
            650,
            655,
            660,
            665,
            670,
            675,
            680,
            685,
            690,
            695,
            700,
            705,
            710,
            999999,
        ],

        # collection efficiency (with factor 0.9) correction
        0.01 *
        0.9 *
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.1,
            0.1,
            0.9,
            1.8,
            4.6,
            7.4,
            11.0,
            15.0,
            18.0,
            21.0,
            22.0,
            23.0,
            24.0,
            25.0,
            25.0,
            26.0,
            26.0,
            26.0,
            26.0,
            27.0,
            26.0,
            26.0,
            26.0,
            26.0,
            26.0,
            25.0,
            25.0,
            25.0,
            25.0,
            24.0,
            24.0,
            24.0,
            23.0,
            22.0,
            22.0,
            21.0,
            20.0,
            19.0,
            19.0,
            19.0,
            18.0,
            18.0,
            18.0,
            17.0,
            16.0,
            15.0,
            14.0,
            12.0,
            11.0,
            10.0,
            9.5,
            8.7,
            8.1,
            7.6,
            7.2,
            6.8,
            6.4,
            6.0,
            5.6,
            5.1,
            4.8,
            4.4,
            4.0,
            3.6,
            3.3,
            3.0,
            2.7,
            2.3,
            2.1,
            1.9,
            1.7,
            1.5,
            1.4,
            1.2,
            1.0,
            0.9,
            0.7,
            0.6,
            0.5,
            0.4,
            0.3,
            0.2,
            0.1,
            0.0,
            0.0,
        ],
        NoBoundaries(),
    ),

    # using eps(), max angle and two zero values in the second vector to
    # get exact zeros outside of the upper boundary defined in Jpp (i.e. >=0.4)
    LinearInterpolator(
        [
            -1.00,
            -0.95,
            -0.90,
            -0.85,
            -0.80,
            -0.75,
            -0.70,
            -0.65,
            -0.60,
            -0.55,
            -0.50,
            -0.45,
            -0.40,
            -0.35,
            -0.30,
            -0.25,
            -0.20,
            -0.15,
            -0.10,
            -0.05,
            0.00,
            0.05,
            0.10,
            0.15,
            0.20,
            0.25,
            0.30,
            0.35,
            0.40,
            0.40 + eps(),
            1.0,
        ],
        [
            1.621,
            1.346,
            1.193,
            1.073,
            0.973,
            0.877,
            0.790,
            0.711,
            0.640,
            0.575,
            0.517,
            0.450,
            0.396,
            0.341,
            0.295,
            0.249,
            0.207,
            0.166,
            0.128,
            0.095,
            0.065,
            0.038,
            0.017,
            0.006,
            0.003,
            0.002,
            0.001,
            0.001,
            0.001,
            0.000,
            0.000,
        ],
        NoBoundaries(),
    ),
)