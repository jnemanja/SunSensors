abstract type AngleDeterminationMethod end

"""
    SpatialMoment

Spatial moment method for determining spot position.
"""
struct SpatialMoment <: AngleDeterminationMethod end

"""
    spot_position(readings, sensor, :SpatialMoment)

Calculates a light spot position on the Profile sensors sensing surface,
using spatial moment method.
"""
function spot_position(
    signal::Array{T}, sensor::ProfileSensor, ::SpatialMoment
) where {T<:Number}
    sensor_length = size(sensor)

    shift = []
    for axis in 1:2
        shift_px =
            sum(signal[:, axis] .* collect(0:((sensor.resolution[axis]) - 1))) /
            sum(signal[:, axis])
        push!(shift, sensor_length[axis] * (shift_px / (sensor.resolution[axis] - 1) - 0.5))
    end

    return shift
end

"""
    QuadraticFitting(window_size)

Quadratic fitting method with configurable window size for determining spot position.
"""
struct QuadraticFitting <: AngleDeterminationMethod
    window::UInt8
end

"""
    spot_position(readings, sensor, :QuadraticFitting)

Calculates a light spot position on the Profile sensors sensing surface,
using quadratic fitting method.
"""
function spot_position(
    signal::Array{T}, sensor::ProfileSensor, method::QuadraticFitting
) where {T<:Number}
    ## Least squares fitting of the quadratic function y(x) = a + bx + cx^2.
    @assert method.window % 2 == 1  ## Works only for odd sized windows.
    n = Integer((method.window + 1) / 2)
    N = Integer((method.window - 1) / 2)
    sensor_length = size(sensor)

    shift = []
    for axis in 1:2
        k = findmax(signal[:, axis])[2]
        if k < N + 1
            k = N + 1
        end
        y = signal[:, axis][(k - N):(k + N)]

        xs1 = (1:(method.window)) .- n
        xs0 = xs1 .^ 0
        xs2 = xs1 .^ 2
        y0 = sum(y .* xs0)
        y1 = sum(y .* xs1)
        y2 = sum(y .* xs2)
        a0 = sum(xs0 .* xs0)
        a2 = sum(xs0 .* xs2)
        b1 = sum(xs1 .* xs1)
        c0 = sum(xs2 .* xs0)
        c2 = sum(xs2 .* xs2)

        c = (y0 * a2 - y2 * a0) / (c0 * a2 - c2 * a0)
        b = y1 / b1
        x = -b / 2c

        shift_px = k - 1 + x
        push!(shift, sensor_length[axis] * (shift_px / (sensor.resolution[axis] - 1) - 0.5))
    end

    return shift
end

"""
    spot_position(signal, sensor)

Default spot position calculation for profile sensor using spatial moment method.
"""
spot_position(signal, sensor::ProfileSensor) =
    spot_position(signal, sensor, SpatialMoment())