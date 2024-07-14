# Tests Hamamatsu sensor models for S5990 and S5991.
# The test measures and compares dark current, resolution and position detection errors reported in
# the sensors datasheet:
# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s5990-01_etc_kpsd1010e.pdf

using Statistics
using SunSensors
using Plots

include("../scripts/hamamatsu_s5990.jl")
include("../scripts/hamamatsu_s5991.jl")

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
    x = ((s[2] + s[3]) - (s[1] + s[4])) / (2 * sum(s)) * L
    y = ((s[1] + s[2]) - (s[3] + s[4])) / (2 * sum(s)) * L
    return (x, y)
end

## Dark current
currents_s5990 = []
currents_s5991 = []
for _ in 1:1000
    for t in 1:5
        exposure_time = t * 1.0e-3u"s"
        push!(
            currents_s5990,
            sum(SunSensors.read(s5990, NoSignal(), exposure_time) / exposure_time),
        )
        push!(
            currents_s5991,
            sum(SunSensors.read(s5991, NoSignal(), exposure_time) / exposure_time),
        )
    end
end
dark_current_s5990 = sum(currents_s5990) / length(currents_s5990)
dark_current_s5991 = sum(currents_s5991) / length(currents_s5991)
println("Dark currents from datasheet S5990 0.5 nA, S5991 1 nA")
println("Dark current S5990 $dark_current_s5990, S5991 $dark_current_s5991")

## Resolution
N = 5000
t = 0.5e-3u"s"
L_s5990 = SunSensors.size(s5990)
L_s5991 = SunSensors.size(s5991)

# Photocurrent should be 1uA.
ls_photocurrent_1uA = RayAngles(0.0u"°", 0.0u"°", 75.0u"W/m^2")
signal = OpticalSignal(lightspot_mask(0.0, 0.0), ls_photocurrent_1uA)
readouts_s5990 = [(SunSensors.read(s5990, signal, t)) / t for _ in 1:N]
readouts_s5991 = [(SunSensors.read(s5991, signal, t)) / t for _ in 1:N]
photocurrent_s5990 = sum(sum.(readouts_s5990)) / N
photocurrent_s5991 = sum(sum.(readouts_s5991)) / N
println("Photocurrents used in datasheet measurements S5990 1 uA, S5991 1 uA")
println("Photocurrents measured $photocurrent_s5990 $photocurrent_s5991")

# Calculate expected resolution from datasheet information, by removing the circuit noise.
circuit_input_noise = 1.0e-6u"V"
measurement_frequency = 1.0e3u"Hz"
measurement_photocurrent = 1.0e-6u"A"
equivalent_input_voltage_noise = circuit_input_noise / sqrt(measurement_frequency)

equivalent_input_current_noise_s5990 =
    equivalent_input_voltage_noise / s5990.photodetector.resistance *
    sqrt(measurement_frequency)
measured_position_resolution_s5990 = 0.7e-6u"m"
measured_noise_s5990 =
    measured_position_resolution_s5990 / L_s5990[1] * measurement_photocurrent
expected_noise_s5990 = sqrt(measured_noise_s5990^2 - equivalent_input_current_noise_s5990^2)
resolution_s5990 = L_s5990[1] * expected_noise_s5990 / measurement_photocurrent

equivalent_input_current_noise_s5991 =
    equivalent_input_voltage_noise / s5991.photodetector.resistance *
    sqrt(measurement_frequency)
measured_position_resolution_s5991 = 1.5e-6u"m"
measured_noise_s5991 =
    measured_position_resolution_s5991 / L_s5991[1] * measurement_photocurrent
expected_noise_s5991 = sqrt(measured_noise_s5991^2 - equivalent_input_current_noise_s5991^2)
resolution_s5991 = L_s5991[1] * expected_noise_s5991 / measurement_photocurrent
println("Expected resolution S5990 $resolution_s5990 S5991 $resolution_s5991")

function measure(x, y, sensor)
    signal = OpticalSignal(lightspot_mask(x, y), ls_photocurrent_1uA)
    out = SunSensors.read(sensor, signal, t)
    pos = spot_position(out, sensor)
    return pos
end

measure_s5990(x, y) = measure(x, y, s5990)
measure_s5991(x, y) = measure(x, y, s5991)

steps = 20
samples = 20

positions_s5990 = []
accuracies_s5990 = []
precisions_s5990 = []
positions_s5991 = []
accuracies_s5991 = []
precisions_s5991 = []

for i in 0:(steps - 1)
    for j in 0:(steps - 1)
        pos_x = (-0.5 + i / (steps - 1)) * ustrip(L_s5990[1]) * 0.9
        pos_y = (-0.5 + j / (steps - 1)) * ustrip(L_s5990[2]) * 0.9
        a, p = measure_accuracy_precision(measure_s5990, (pos_x, pos_y), 0.05e-3, samples)
        push!(positions_s5990, (pos_x, pos_y))
        push!(accuracies_s5990, a)
        push!(precisions_s5990, p)
        pos_x = (-0.5 + i / (steps - 1)) * ustrip(L_s5991[1]) * 0.9
        pos_y = (-0.5 + j / (steps - 1)) * ustrip(L_s5991[2]) * 0.9
        a, p = measure_accuracy_precision(measure_s5991, (pos_x, pos_y), 0.05e-3, samples)
        push!(positions_s5991, (pos_x, pos_y))
        push!(accuracies_s5991, a)
        push!(precisions_s5991, p)
    end
end

println(
    "Maximum and typical position detection errors in datasheet for S5990 150 μm and 70 μm"
)
println(
    "Maximum and typical position detection errors in datasheet for S5991 250 μm and 150 μm"
)

avg_precision_s5990 = mean(precisions_s5990)
avg_precision_s5991 = mean(precisions_s5991)
println(
    "Average measured resolutions S5990 $avg_precision_s5990 S5991 $avg_precision_s5991"
)
avg_accuracy_s5990 = sqrt(
    sum(mean.([first.(accuracies_s5990), last.(accuracies_s5990)]) .^ 2)
)
avg_accuracy_s5991 = sqrt(
    sum(mean.([first.(accuracies_s5991), last.(accuracies_s5991)]) .^ 2)
)
println(
    "Average measured detection errors S5990 $avg_accuracy_s5990 S5991 $avg_accuracy_s5991"
)
t_accuracy_s5990 = sqrt.(sum(avg_accuracy_s5990 .^ 2))
t_accuracy_s5991 = sqrt.(sum(avg_accuracy_s5991 .^ 2))
println(
    "Average measured total detection errors S5990 $t_accuracy_s5990 S5991 $t_accuracy_s5991",
)

pos_s5990 = [(-0.5 + i / (steps - 1)) * ustrip(L_s5990[1]) * 0.9 for i in 0:(steps - 1)]
pos_s5991 = [(-0.5 + i / (steps - 1)) * ustrip(L_s5991[1]) * 0.9 for i in 0:(steps - 1)]
output_dir = "data/output/test/"
heatmap(
    pos_s5990,
    pos_s5990,
    reshape(precisions_s5990 .* 1e6, steps, steps);
    label="S5990 resolution",
    xlabel="position [m]",
    ylabel="position [m]",
    colorbar_title="precision [μm]",
)
savefig(output_dir * "s5990_precision.svg")
heatmap(
    pos_s5991,
    pos_s5991,
    reshape(precisions_s5991 * 1e6, steps, steps);
    label="S5991 resolution",
    xlabel="position [m]",
    ylabel="position [m]",
    colorbar_title="precision [μm]",
)
savefig(output_dir * "s5991_precision.svg")
heatmap(
    pos_s5990,
    pos_s5990,
    reshape(sqrt.(sum.(broadcast(x -> x .^ 2, accuracies_s5990))) * 1e6, steps, steps);
    label="S5990 position detection error",
    xlabel="position [m]",
    ylabel="position [m]",
    colorbar_title="accuracy [μm]",
)
savefig(output_dir * "s5990_accuracy.svg")
heatmap(
    pos_s5991,
    pos_s5991,
    reshape(sqrt.(sum.(broadcast(x -> x .^ 2, accuracies_s5991))) * 1e6, steps, steps);
    label="S5991 position detection error",
    xlabel="position [m]",
    ylabel="position [m]",
    colorbar_title="accuracy [μm]",
)
savefig(output_dir * "s5991_accuracy.svg")

println("Accuracy and precision plots for S5990 and S5991 models saved in $output_dir")
