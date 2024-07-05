# Tests Hamamatsu sensor models for S3931 and S3932.
# The test measures and compares dark current, resolution and position detection errors reported in
# the sensors datasheet:
# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s3931_etc_kpsd1002e.pdf

using SunSensors
using Plots

include("../scripts/hamamatsu_s3931.jl")
include("../scripts/hamamatsu_s3932.jl")

# Creates the mask with the pinhole at the (x, y) position.
function lightspot_mask(x, y)
    return PositionedMask(
        BasicMask(
            [PinholeFeature(100.0e-6u"m")], [Placement(1, Point(x, y), 0.0)], 1.0e-8u"m"
        ),
        (0, 0, 1.0e-3),
        RotZ(0),
    )
end

function get_pos(signal, L)
    s = ustrip.(signal)
    return (s[2] - s[1]) / sum(s) / 2 * L
end

## Dark current
currents_s3931 = []
currents_s3932 = []
for _ in 1:1000
    for t in 1:5
        exposure_time = t * 1.0e-3u"s"
        push!(
            currents_s3931,
            sum(SunSensors.read(s3931, NoSignal(), exposure_time) / exposure_time),
        )
        push!(
            currents_s3932,
            sum(SunSensors.read(s3932, NoSignal(), exposure_time) / exposure_time),
        )
    end
end
dark_current_s3931 = uconvert(u"nA", sum(currents_s3931) / length(currents_s3931))
dark_current_s3932 = uconvert(u"nA", sum(currents_s3932) / length(currents_s3932))
println("Dark currents from datasheet S3931 0.15 nA, S3932 0.2 nA")
println("Dark currents measured S3931 $dark_current_s3931, S3932 $dark_current_s3932")

## Resolution
ls_photocurrent_1uA = RayAngles(0.0u"°", 0.0u"°", 102.0u"W/m^2")
signal = OpticalSignal(lightspot_mask(0.0, 0.0), ls_photocurrent_1uA)

# Photocurrent should be 1uA.
N = 5000
t = 0.5e-3u"s"
L_s3931 = 6.0e-3u"m"
L_s3932 = 12.0e-3u"m"

# Calculate expected resolution from datasheet information, by removing the circuit noise.
circuit_input_noise = 1.0e-6u"V"
measurement_frequency = 1.0e3u"Hz"
measurement_photocurrent = 1.0e-6u"A"
equivalent_input_voltage_noise = circuit_input_noise / sqrt(measurement_frequency)

equivalent_input_current_noise_s3931 =
    equivalent_input_voltage_noise / s3931.photodetector.resistance *
    sqrt(measurement_frequency)
measured_position_resolution_s3931 = 0.2e-6u"m"
measured_noise_s3931 =
    measured_position_resolution_s3931 / L_s3931 * measurement_photocurrent
expected_noise_s3931 = sqrt(measured_noise_s3931^2 - equivalent_input_current_noise_s3931^2)
resolution_s3931 = L_s3931 * expected_noise_s3931 / measurement_photocurrent

equivalent_input_current_noise_s3932 =
    equivalent_input_voltage_noise / s3932.photodetector.resistance *
    sqrt(measurement_frequency)
measured_position_resolution_s3932 = 0.3e-6u"m"
measured_noise_s3932 =
    measured_position_resolution_s3932 / L_s3932 * measurement_photocurrent
expected_noise_s3932 = sqrt(measured_noise_s3932^2 - equivalent_input_current_noise_s3932^2)
resolution_s3932 = L_s3932 * expected_noise_s3932 / measurement_photocurrent
println("Expected resolutions S3931 $resolution_s3931 S3932 $resolution_s3932")

readouts_s3931 = [(SunSensors.read(s3931, signal, t)) / t for _ in 1:N]
readouts_s3932 = [(SunSensors.read(s3932, signal, t)) / t for _ in 1:N]

# Resolution is approximately the standard deviation.
positions_s3931 = [(s[2] - s[1]) / (2 * sum(s)) * L_s3931 for s in readouts_s3931]
mean_s3931 = sum(positions_s3931) / N
sd_s3931 = sqrt(sum([(p - mean_s3931)^2 for p in positions_s3931]) / N)
println("S3931 mean $mean_s3931, std $sd_s3931")

positions_s3932 = [(s[2] - s[1]) / (2 * sum(s)) * L_s3932 for s in readouts_s3932]
mean_s3932 = sum(positions_s3932) / N
sd_s3932 = sqrt(sum([(p - mean_s3932)^2 for p in positions_s3932]) / N)
println("S3932 mean $mean_s3932, std $sd_s3932")

photocurrent_s3931 = sum(sum.(readouts_s3931)) / N
photocurrent_s3932 = sum(sum.(readouts_s3932)) / N
println("Photocurrents used in datasheet measurements S3931 1 uA, S3932 1 uA")
println("Photocurrents measured S3931 $photocurrent_s3931, S3932 $photocurrent_s3932")

## Position detection error
readouts_s3931 = []
readouts_s3932 = []
steps = 100
ls = RayAngles(0.0u"°", 0.0u"°")
for step in 1:steps
    k = -0.5 + step / steps
    x_s3931 = k * L_s3931
    s_s3931 = SunSensors.read(s3931, OpticalSignal(lightspot_mask(ustrip(x_s3931), 0.0), ls), t) / t
    push!(readouts_s3931, (x_s3931, get_pos(s_s3931, L_s3931)))
    x_s3932 = k * L_s3932
    s_s3932 = SunSensors.read(s3932, OpticalSignal(lightspot_mask(ustrip(x_s3932), 0.0), ls), t) / t
    push!(readouts_s3932, (x_s3932, get_pos(s_s3932, L_s3932)))
end

p3931 = plot(
    first.(readouts_s3931),
    last.(readouts_s3931) .- first.(readouts_s3931);
    ylims=(-50e-6, 50e-6),
)

p3932 = plot(
    first.(readouts_s3932),
    last.(readouts_s3932) .- first.(readouts_s3932);
    ylims=(-50e-6, 50e-6),
)

# Plot to be compared with the datasheet figures on page 4, "Examples of position detectability".
plot(p3931, p3932; layout=(2, 1))
