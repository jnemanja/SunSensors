"""
    OpticalSignal

Models the incident light on the sensing surface of the photodetector.
"""
abstract type OpticalSignal end

"""
    NoSignal

Models absence of optical signal on the sensing  surface of the photodetector.
"""
struct NoSignal <: OpticalSignal end

"""
    PatternOpticalSignal(pattern, power)

Models the incident light pattern on the sensing surface of the photodetector.
"""
struct PatternOpticalSignal <: OpticalSignal
    # The pattern of the light.
    pattern::GeometrySet

    # The power of the light.
    power_density::typeof(1.0u"W/m^2")
end

function extend_point(p::Point)
    p1, p2 = (to(p)...,)
    p3 = zero(p1)
    return Point(p1, p2, p3)
end

function projectray(point::Point, ray::Point)
    vp = to(point)
    vr = to(ray)
    p = Point((vp[1:2] .- (vp[3] / vr[3]) * vr[1:2])...)
    return p
end

"""
    OpticalSignal(mask, light_source)

Constructs the OpticalSignal using a light source and a mask.
"""
function OpticalSignal(mask::PositionedMask, ls::LightSource)
    # Construct a pattern of the projected light throught the mask.
    rot_m = Rotate(mask.rotation)
    ray_vector = Point(direction_vector(ls)...)
    features = []
    for placement in mask.mask.feature_placements
        rot_f = Rotate(placement.rotation)
        feature = Translate(to(placement.position)...)(
            rot_f(geometry(mask.mask.feature_set[placement.feature_index]))
        )
        feature_3d = Ngon([extend_point(p) for p in pointify(feature)]...)

        t = mask.mask.thickness / 2
        feature_bottom = Translate(mask.shift...)(rot_m(Translate(0.0u"m", 0.0u"m", -t)(feature_3d)))
        feature_top = Translate(mask.shift...)(rot_m(Translate(0.0u"m", 0.0u"m", t)(feature_3d)))
        feautre_projection_bottom = Ngon(
            (projectray(p, ray_vector) for p in pointify(feature_bottom))...
        )
        feautre_projection_top = Ngon(
            (projectray(p, ray_vector) for p in pointify(feature_top))...
        )
        new = intersection(feautre_projection_top, feautre_projection_bottom).geom
        push!(features, new)
    end
    filter!((!isnothing), features)

    # Calculate the power density of the optical signal.
    if length(features) > 0
        ray_vector = direction_vector(ls)
        light_power_nominal = power(ls)
        light_power = ray_vector[3] * light_power_nominal

        return PatternOpticalSignal(GeometrySet(features), light_power)
    else
        return NoSignal()
    end
end
