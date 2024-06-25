"""
    Placement(feature_index, position, rotation)

Placement of a feature on the mask, characterized by the feature's index in the feature set, its
position on the mask and the its rotation.
"""
struct Placement
    # Index of the feature in the feature set.
    feature_index::UInt8

    # Position of the feature on the mask.
    position::Point

    # Rotation of the feature.
    rotation::typeof(1.0u"°")
end

"""
    Feature

Feature on the mask.
"""
abstract type Feature end

"""
    geometry(feature)

Returns a geometry of the feature.
"""
function geometry end

"""
    PinholeFeautre(radius, sides)

Pihnole feautre with circular geometry (polygon with given sides) and given radius.
"""
struct PinholeFeature <: Feature
    # Radius of the pinhole.
    radius::typeof(1.0u"m")

    # Number of sides of the representing polygon.
    sides::UInt32
end

"""
    PinholeFeature(radius)

Constructs the PinholeFeature of the given radius with a 20 sided polygon.
"""
PinholeFeature(radius::typeof(1.0u"m")) = PinholeFeature(radius, 20)

"""
    geometry(pinhole_feature)

Returns geometry of the PinholeFeature.
"""
function geometry(f::PinholeFeature)
    angle_step = 2 * pi / f.sides
    return Ngon{20}(
        [(sin(i * angle_step), cos(i * angle_step)) .* f.radius for i in 1:(f.sides)]...
    )
end

"""
    LFeature(width, height, thickness)

L-shaped feature with the given width, height and thickness.
"""
struct LFeature <: Feature
    # Width of the L shape.
    width::typeof(1.0u"m")

    # Height of the L shape.
    height::typeof(1.0u"m")

    # Thickness of the L shape.
    thickness::typeof(1.0u"m")
end

"""
    geometry(l_feature)

Returns geometry of the LFeature.
"""
function geometry(f::LFeature)
    return Translate(-f.width / 2, -f.height / 2)(
        Ngon(
            (0.0, 0.0),
            (f.width, 0.0),
            (f.width, f.thickness),
            (f.width - f.thickness, f.thickness),
            (f.width - f.thickness, f.height),
            (0.0, f.height),
        ),
    )
end

"""
    RectangleFeature(width, height)

Rectangular feature of the given width and height.
"""
struct RectangleFeature <: Feature
    width::typeof(1.0u"m")
    height::typeof(1.0u"m")
end

"""
    geometry(l_feature)

Returns geometry of the RectangleFeature.
"""
function geometry(f::RectangleFeature)
    x = f.width / 2
    y = f.height / 2
    return Quadrangle((-x, -y), (x, -y), (x, y), (-x, y))
end
