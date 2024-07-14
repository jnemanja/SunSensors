# Tests Hamamatsu sensor models for S3931 and S3932.
# The test measures and compares dark current, resolution and position detection errors reported in
# the sensors datasheet:
# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s3931_etc_kpsd1002e.pdf

using Statistics
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
N = 5000
t = 0.5e-3u"s"
L_s3931 = SunSensors.size(s3931)[1]
L_s3932 = SunSensors.size(s3932)[1]

# Photocurrent should be 1uA.
ls_photocurrent_1uA = RayAngles(0.0u"°", 0.0u"°", 102.0u"W/m^2")
signal = OpticalSignal(lightspot_mask(0.0, 0.0), ls_photocurrent_1uA)
readouts_s3931 = [(SunSensors.read(s3931, signal, t)) / t for _ in 1:N]
readouts_s3932 = [(SunSensors.read(s3932, signal, t)) / t for _ in 1:N]
photocurrent_s3931 = sum(sum.(readouts_s3931)) / N
photocurrent_s3932 = sum(sum.(readouts_s3932)) / N
println("Photocurrents used in datasheet measurements S3931 1 uA, S3932 1 uA")
println("Photocurrents measured S3931 $photocurrent_s3931, S3932 $photocurrent_s3932")

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

function measure(x, y, sensor)
    signal = OpticalSignal(lightspot_mask(x, y), ls_photocurrent_1uA)
    out = SunSensors.read(sensor, signal, t)
    pos = spot_position(out, sensor)
    return (pos, y)
end

measure_s3931(x, y) = measure(x, y, s3931)
measure_s3932(x, y) = measure(x, y, s3932)

steps = 100
samples = 50

positions_s3931 = []
accuracies_s3931 = []
precisions_s3931 = []
positions_s3932 = []
accuracies_s3932 = []
precisions_s3932 = []

for i in 0:(steps - 1)
    pos = (-0.5 + i / (steps - 1)) * ustrip(L_s3931) * 0.9
    a, p = measure_accuracy_precision(measure_s3931, (pos, 0.0), 0.05e-3, samples)
    push!(positions_s3931, pos)
    push!(accuracies_s3931, a)
    push!(precisions_s3931, p)
    pos = (-0.5 + i / (steps - 1)) * ustrip(L_s3932) * 0.9
    a, p = measure_accuracy_precision(measure_s3932, (pos, 0.0), 0.05e-3, samples)
    push!(positions_s3932, pos)
    push!(accuracies_s3932, a)
    push!(precisions_s3932, p)
end
avg_precision_s3931 = mean(precisions_s3931)
avg_precision_s3932 = mean(precisions_s3932)
println(
    "Average measured resolutions S3931 $avg_precision_s3931 S3932 $avg_precision_s3932"
)
avg_accuracy_s3931 = mean(first.(accuracies_s3931))
avg_accuracy_s3932 = mean(first.(accuracies_s3932))
println(
    "Average measured detection errors S3931 $avg_accuracy_s3931 S3932 $avg_accuracy_s3932"
)

output_dir = "data/output/test/"
plot(
    positions_s3931,
    precisions_s3931;
    label="S3931 resolution",
    xlabel="position [m]",
    ylabel="precision [m]",
)
savefig(output_dir * "s3931_precision.svg")
plot(
    positions_s3932,
    precisions_s3932;
    label="S3932 resolution",
    xlabel="position [m]",
    ylabel="precision [m]",
)
savefig(output_dir * "s3932_precision.svg")
plot(
    positions_s3931,
    first.(accuracies_s3931);
    label="S3931 position detection error",
    xlabel="position [m]",
    ylabel="accuracy [m]",
)
savefig(output_dir * "s3931_accuracy.svg")
plot(
    positions_s3932,
    first.(accuracies_s3932);
    label="S3932 position detection error",
    xlabel="position [m]",
    ylabel="accuracy [m]",
)
savefig(output_dir * "s3932_accuracy.svg")

println("Accuracy and precision plots for S3931 and S3932 models saved in $output_dir")
