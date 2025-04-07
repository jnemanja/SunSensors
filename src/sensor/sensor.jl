"""
    Sensor

An assembly of photodetectors and a signal converter. Sensors provide reading output of the
incident optical signal with read(sensor, optical_signal, exposure_time) function.
"""
abstract type Sensor end

"""
    read(sensor, optical_signal, exposure_time)

Reads sensor output with an incident optical signal, with configurable exposure time.
"""
function read end

"""
    CcdSensor(resolution, pixel)

CCD sensor model with configurable resolution of two-dimensional photodetector array.
Pixels of the sensor are characterized by the same photodetector model.
"""
struct CcdSensor <: Sensor
    # Resolution of the CCD sensor.
    resolution::Tuple{Number,Number}

    # Photodetector characterizing sensor's pixels.
    pixel::Photodetector
end

# Helper function for identifying the pixel indexes in the 2D array given a coordinate.
# Clamps the result to the closest edge if the coordinate is outside the sensors dimensions.
function pixelindex(position, sensorsize, pixelsize, resolution)
    i = Integer.(floor.((2 .* position .+ sensorsize) ./ (2 .* pixelsize))) .+ (1, 1)
    return clamp.(i, range.(1, Signed.(resolution)))
end

# Holds references to features which contribute to signal received by a single pixel.
struct PixelFeatures{N}
    count::Int64
    features::SVector{N,UInt16}
end

function PixelFeatures{N}(feature) where {N}
    return PixelFeatures{N}(1, SVector{N}(feature, repeat([0], N - 1)...))
end

function zero(::Type{PixelFeatures{N}}) where {N}
    return PixelFeatures{N}(0, SVector{N}(repeat([0], N)...))
end

function combine(pf1::PixelFeatures{N}, pf2::PixelFeatures{N}) where {N}
    return PixelFeatures{N}(
        pf1.count + pf2.count,
        vcat(pf1.features[1:(pf1.count)], pf2.features[1:(N - pf1.count)]),
    )
end

# read function for CCD sensors with parametrized size for PixelFeature.
# For optimizing on the expected number of features per pixel.
function read(
    sensor::CcdSensor, os::OpticalSignal, exposure_time::typeof(1.0u"s"), N::Int64
)
    pixel_features = zeros(PixelFeatures{N}, sensor.resolution...)

    pixelbox = Meshes.boundingbox(sensor.pixel.geometry)
    pixelsize = pixelbox.max - pixelbox.min
    sensor_size = pixelsize .* sensor.resolution
    sensor_bottomleft = sensor_size ./ -2
    sensor_topright = sensor_size ./ 2
    sensor_box = Box(Point(sensor_bottomleft...), Point(sensor_topright...))

    if typeof(os) !== NoSignal
        # For each feature placement, find the pixels it might affect.
        @inbounds for (pi, p) in enumerate(os.pattern)
            bbox = Meshes.boundingbox(p)

            if !intersects(p, sensor_box)
                continue
            end

            # Limit the bounding box to the sensor surface.
            cornerbottomleft = clamp.(to(bbox.min), sensor_bottomleft, sensor_topright)
            cornertopright = clamp.(to(bbox.max), sensor_bottomleft, sensor_topright)

            pixelbottomleft =
                cornerbottomleft .-
                (pixelsize .+ (sensor_size ./ 2 .+ cornerbottomleft) .% pixelsize)
            pixeltopright =
                cornertopright .+
                (pixelsize .- (sensor_size ./ 2 .+ cornertopright) .% pixelsize)
            pixelbottomleft_index = pixelindex(
                pixelbottomleft, sensor_size, pixelsize, sensor.resolution
            )
            pixeltopright_index = pixelindex(
                pixeltopright .- pixelsize ./ 2, sensor_size, pixelsize, sensor.resolution
            )
            x_range = pixelbottomleft_index[1]:pixeltopright_index[1]
            y_range = pixelbottomleft_index[2]:pixeltopright_index[2]
            pixel_features[x_range, y_range] .=
                combine.(
                    pixel_features[x_range, y_range],
                    repeat([PixelFeatures{N}(pi)], length(x_range), length(y_range)),
                )
        end
    end

    out = zeros(typeof(1.0u"C"), sensor.resolution[2], sensor.resolution[1])
    # For each pixel, build a mask with only relevant feature placements,
    # and read that pixel.
    @inbounds for y in 1:sensor.resolution[2]
        @inbounds for x in 1:sensor.resolution[1]
            yy = (y - 0.5 - sensor.resolution[2] / 2) * pixelsize[2]
            xx = (x - 0.5 - sensor.resolution[1] / 2) * pixelsize[1]

            fs = []
            for i in 1:(pixel_features[x, y].count)
                push!(
                    fs,
                    Translate(-xx, -yy)(os.pattern.geoms[pixel_features[x, y].features[i]]),
                )
            end
            if length(fs) == 0
                new_os = NoSignal()
            else
                new_os = PatternOpticalSignal(GeometrySet(fs), os.power_density)
            end
            out[y, x] = charge(sensor.pixel, new_os, exposure_time)
        end
    end

    return out
end

"""
    read(ccd_sensor, optical_signal, exposure_time)

Reads the optical signal with the CCD sensor, after the provided exposure time.
"""
function read(s::CcdSensor, os::OpticalSignal, t::typeof(1.0u"s"))
    # Assume the CCD sensor pixels are sufficiently smaller than the optical sgianl featurs and
    # four are to be read by a pixel at maximum.
    return read(s, os, t, 4)
end

function size(s::CcdSensor)
    bbox = boundingbox(s.pixel.geometry)
    return (bbox.max - bbox.min) .* s.resolution
end

"""
    ProfileSensor(resolution, horizontal_pixel, vertical_pixel)

Profile sensor model with configurable resolution and pixels for horizontal and vertical arrays.
Consists of two orthogonaly oriented one-dimensional photodetector arrays overlayed over the same
sensing surface.

This constructor allows unrealistic sensor geometries if the resolution and pixel geometries do not
form a proper rectangle.
"""
struct ProfileSensor <: Sensor
    # Resolution of the profile sensor.
    resolution::Tuple{Unsigned,Unsigned}

    # Photodetector characterizing sensors's horizontal pixels.
    horizontal_pixel::Photodetector

    # Photodetector characterizing sensors's vertical pixels.
    vertical_pixel::Photodetector
end

"""
    ProfileSensor(resolution, horizontal_pixel)

Constructs profile sensor model with provided resolution and horizontal pixel.
The vertical pixel is produced by rotating the horizontal pixel by 90 degrees.
"""
function ProfileSensor(resolution::Tuple{Number,Number}, horizontal_pixel::Photodetector)
    # Construct vertical pixel by rotating the horizontal.
    if typeof(horizontal_pixel) <: IdealPhotodetector
        vertical_pixel = IdealPhotodetector(
            Rotate(pi / 2)(horizontal_pixel.geometry),
            horizontal_pixel.responsivity,
            horizontal_pixel.dark_power_density,
        )
    elseif typeof(horizontal_pixel) <: NoisyPhotodetector
        vertical_pixel = NoisyPhotodetector(
            Rotate(pi / 2)(horizontal_pixel.geometry),
            horizontal_pixel.responsivity,
            horizontal_pixel.dark_power_density,
            horizontal_pixel.temperature,
            horizontal_pixel.resistance,
        )
    else
        error(
            "Unsupported photodetector type.
            Consider manually constructing photodetector for vertical pixel and using
            ProfileSensor(resolution, horizontal_pixel, vertical_pixel, adc) constructor."
        )
    end

    return ProfileSensor(resolution, horizontal_pixel, vertical_pixel)
end

"""
    read(profile_sensor, optical_signal, exposure_time)

Reads the optical signal with the profile sensor, after the provided exposure time.
"""
function read(sensor::ProfileSensor, os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    ccd_horizontal = CcdSensor((sensor.resolution[1], 1), sensor.horizontal_pixel)
    ccd_vertical = CcdSensor((1, sensor.resolution[2]), sensor.vertical_pixel)
    # Pixels of the profile sensors are elongated and possible to sense large number of features.
    out_horizontal = read(ccd_horizontal, os, exposure_time, 16)
    out_vertical = read(ccd_vertical, os, exposure_time, 16)

    return hcat(out_horizontal', out_vertical)
end

function size(s::ProfileSensor)
    v_bbox = boundingbox(s.vertical_pixel.geometry)
    h_bbox = boundingbox(s.horizontal_pixel.geometry)
    w = max((h_bbox.max - h_bbox.min)...)
    h = max((v_bbox.max - v_bbox.min)...)
    return (w, h)
end

"""
    PsdSensor(photodetector, transfer_coefficinets)

Models PSD (Position Sensitive Device) sensor with configurable photodetector, signal_converter and
transfer coefficinets. Transfer coefficinets is a matrix with a row for each of the sensor outputs,
and 5 columns for (I) constant offsets, (II and III) linear and (IV and V) squared coordinates
coefficinets.
"""
struct PsdSensor <: Sensor
    # Photodetector for the PDS sensor's sensing surface.
    photodetector::Photodetector

    # Transfer coefficinets matrix with a row for each output value and 5 columns.
    # First column contains constant offsets, second and third column are linear factors for X and
    # Y coordinate, and fourth and fifth are coefficinets for squared X and Y coordinates.
    transfer_coeffs::Matrix{Float64}
end

"""
    read(psd_sensor, optical_signal, exposure_time)

Reads the optical signal with the PSD sensor, after the provided exposure time.
"""
function read(sensor::PsdSensor, os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    c = charge(sensor.photodetector, os, exposure_time)

    position = centroid(sensor.photodetector, os, exposure_time)
    if isnothing(position)
        position = Point(0.0, 0.0)
    end

    constant = sensor.transfer_coeffs[:, 1]
    linear = sensor.transfer_coeffs[:, 2:3]
    squared = sensor.transfer_coeffs[:, 4:5]

    p_coords = [uconvert(NoUnits, e / u"m") for e in to(position)]
    outs = constant + linear * p_coords + squared * (p_coords .^ 2 .* sign.(p_coords))
    charges = outs ./ sum(outs) * c

    return charges
end

function size(s::PsdSensor)
    bbox = boundingbox(s.photodetector.geometry)
    return (bbox.max - bbox.min)
end
