module SunSensors

using Reexport: @reexport
using Distributions: Normal
using Distributions: Poisson
using Meshes
@reexport using Meshes: Point
using Rotations: Rotation
@reexport using Rotations: RotX
@reexport using Rotations: RotY
@reexport using Rotations: RotZ
using StaticArrays: SVector
@reexport using Unitful

import Base.zero

include("sensor/light_source.jl")

export RayAngles
export RayVector
export direction_vector
export power

include("sensor/mask_features.jl")

export BasicMask
export PositionedMask

include("sensor/mask.jl")

export PinholeFeature
export Placement
export RectangleFeature

include("sensor/optical_signal.jl")
export NoSignal
export PatternOpticalSignal

include("sensor/photodetector.jl")

export Photodetector
export IdealPhotodetector
export NoisyPhotodetector
export charge
export centroid

include("sensor/signal_converter.jl")

export AAConverter
export ADConverter
export convert

include("sensor/sensor.jl")

export CcdSensor
export ProfileSensor
export PsdSensor
export OpticalSignal
export read

end
