# Tests Hamamatsu sensor model for S9132.
# The test measures and compares dark current and voltage outputs as reported in the sensors datasheet:
# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s9132_kmpd1075e.pdf

using Statistics
using SunSensors
using Plots
using Meshes

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

# Needed to remove errorneous approximations.
import Meshes.atol
Meshes.atol(::Type{Float64}) = 1e-12

include("../scripts/hamamatsu_s9132.jl")

## Dark current
signal = NoSignal()
t = 1e-1u"s"

samples = 5000
pixel_count = 512
readings = [SunSensors.read(s9132, signal, t) for _ in 1:samples]

dark_current = mean(vcat(vcat(readings)...)) / t
bits_8_low = mean(vcat(vcat(adc_8bit_low_gain.(readings))...))
bits_10_low = mean(vcat(vcat(adc_10bit_low_gain.(readings))...))
bits_8_high = mean(vcat(vcat(adc_8bit_high_gain.(readings))...))
bits_10_high = mean(vcat(vcat(adc_10bit_high_gain.(readings))...))

println("Dark current in datasheet 0.2 pA")
println("Dark current $dark_current")
println("Dark output voltages in datasheet are 20 mV for low gain and 100 mV for high gain")
println("Converted to output bits, for low gain 8-bit is 1, 10-bit is 5 and for high gain 8-bit is 7 and 10-bit is 27")
println("Dark bits measured $bits_8_low, $bits_10_low, $bits_8_high, $bits_10_high")
println("Noise significantly biases the outputs, given the signal is clammped above 0")

## Accuracy and precision plots
ls = RayAngles(0.0u"°", 0.0u"°", 75.0u"W/m^2")

function measure(x, y)
    signal = OpticalSignal(lightspot_mask(x, y), ls)
    out = SunSensors.read(s9132, signal, t)
    pos = spot_position(out, s9132)
    return ustrip.(pos)
end

L_s9132 = SunSensors.size(s9132)

steps = 20
samples = 20

positions = []
accuracies = []
precisions = []

for i in 0:(steps - 1)
    for j in 0:(steps - 1)
        pos_x = (-0.5 + i / (steps - 1)) * ustrip(L_s9132[1]) * 0.8
        pos_y = (-0.5 + j / (steps - 1)) * ustrip(L_s9132[2]) * 0.8
        a, p = measure_accuracy_precision(measure, (pos_x, pos_y), 0.05e-3, samples)
        push!(positions, (pos_x, pos_y))
        push!(accuracies, a)
        push!(precisions, p)
    end
end

pos_s9132 = [(-0.5 + i / (steps - 1)) * ustrip(L_s9132[1]) * 0.9 for i in 0:(steps - 1)]
output_dir = "data/output/test/"
heatmap(
    pos_s9132,
    pos_s9132,
    reshape(precisions .* 1e6, steps, steps);
    label="S9132 resolution",
    xlabel="position [m]",
    ylabel="position [m]",
    colorbar_title="precision [μm]",
)
savefig(output_dir * "s9132_precision.svg")

heatmap(
    pos_s9132,
    pos_s9132,
    reshape(sqrt.(sum.(broadcast(x -> x .^ 2, accuracies))) * 1e6, steps, steps);
    label="S9132 position detection error",
    xlabel="position [m]",
    ylabel="position [m]",
    colorbar_title="accuracy [μm]",
)
savefig(output_dir * "s9132_accuracy.svg")
