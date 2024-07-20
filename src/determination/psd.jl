"""
    spot_position(readings, sensor)

Calculates a light spot position on the PSD sensors sensing surface.
Supports one and two dimensional sensors.
"""
function spot_position(readings::Array, sensor::PsdSensor)
    s = sum(readings)
    K = ustrip.(size(sensor)) ./ (2 * s)
    if length(readings) == 2
        if s == 0
            return 0.0
        end
        x = sum(readings .* (-1, 1))
        return x * K[1]
    elseif length(readings) == 4
        if s == 0
            return [0.0, 0.0]
        end
        x = sum(readings .* (-1, 1, 1, -1))
        y = sum(readings .* (1, 1, -1, -1))
        return (x, y) .* K
    end
end
