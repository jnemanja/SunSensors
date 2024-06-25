# Average solar irradiance per unit area at 1 AU.
const solar_power = 1.361e3u"W/m^2"

"""
    LightSource

Models the light source.
"""
abstract type LightSource end

"""
    RayAngles(α, β, power)

Light source defined by the angles and the power.
"""
struct RayAngles <: LightSource
    α::typeof(1.0u"°")
    β::typeof(1.0u"°")
    power::typeof(1.0u"W/m^2")
end

"""
    RayAngles(α, β)

Constructs RayAngles light source using provided angles and the solar power constant.
"""
RayAngles(α::typeof(1.0u"°"), β::typeof(1.0u"°")) = RayAngles(α, β, solar_power)

"""
    direction_vector(ray_angles)

Calculates direction unit vector of the RayAngles light source.
"""
function direction_vector(light_source::RayAngles)
    z = sqrt(1 / (tan(light_source.α)^2 + tan(light_source.β)^2 + 1))
    x = z * tan(light_source.α)
    y = z * tan(light_source.β)
    return SVector(x, y, z)
end

"""
    power(ray_angles)

Power of the RayAngles light source.
"""
power(ls::RayAngles) = ls.power

"""
    RayVector(x, y, z)

Light source defined by the vector with magnitude proportional with the power.
"""
struct RayVector <: LightSource
    x::Float64
    y::Float64
    z::Float64
end

"""
    direction_vector(ray_angles)

Calculates direction unit vector of the RayVector light source.
"""
direction_vector(ls::RayVector) =
    SVector(ls.x, ls.y, ls.z) ./ sqrt((ls.x^2) + (ls.y^2) + (ls.z^2))

"""
    power(ray_angles)

Power of the RayVector light source.
"""
power(ls::RayVector) = sqrt((ls.x^2) + (ls.y^2) + (ls.z^2)) * u"W/m^2"
