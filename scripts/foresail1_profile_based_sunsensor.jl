include("hamamatsu_s9132.jl")

cover_mask_foresail1_profile = PositionedMask(
    BasicMask(
        [PinholeFeature(0.05e-3u"m")], [Placement(1, Point(0.0, 0.0), 0.0)], 1.0e-6u"m"
    ),
    (0.0, 0.0, 2.4e-3),
    RotZ(0.0),
)

function read_profile(ls::LightSource, exposure_time)
    # Sensor is equipped with ND filter with 0.1 transmitance.
    filtered_ls = RayVector((direction_vector(ls) * ustrip(power(ls)) * 0.1)...)
    return s9132_10bit_low_gain(
        OpticalSignal(cover_mask_foresail1_profile, filtered_ls), exposure_time
    )
end
