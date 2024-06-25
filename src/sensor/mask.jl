"""
    Mask

Models the mask for Sun sensors.
"""
abstract type Mask end

"""
    BasicMask(feature_set, feature_placements, thickness)

Mask consisting of a set of features and their placements over the mask's surface with a given
thickness.
"""
struct BasicMask <: Mask
    # A set of features, referred to by the feature placements.
    feature_set::Array{Feature}

    # Feature placements of the features from the set.
    feature_placements::Array{Placement}

    # Thickness of the mask.
    thickness::typeof(1.0u"m")
end

"""
    PositionedMask(mask, shift, rotation)

Models positioning of the mask, using the position shift and rotation.
"""
struct PositionedMask
    # Mask to be positioned.
    mask::Mask

    # Mask's position shift.
    shift::Vec

    # Mask's rotation.
    rotation::Rotation
end

"""
    PositionedMask(mask, height)

Constructs the PositionedMask using a mask and height. Mask is centered in X and Y axes and
parallel with the sensing surface.
"""
function PositionedMask(mask::Mask, height::Number)
    return PositionedMask(mask, Vec(0, 0, height), one(RotMatrix{3}))
end

"""
    PositionedMask(mask, height, rotation)

Constructs the PositionedMask using a mask, height and rotation. Mask is centered in X and Y axes.
"""
function PositionedMask(mask::Mask, height::Number, rotation::Rotation)
    return PositionedMask(mask, Vec(0, 0, height), rotation)
end
