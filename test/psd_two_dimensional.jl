# Tests Hamamatsu sensor models for S5990 and S5991.
# The test measures and compares dark current, resolution and position detection errors reported in
# the sensors datasheet:
# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s5990-01_etc_kpsd1010e.pdf

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
ls_photocurrent_1uA = RayAngles(0.0u"°", 0.0u"°", 75.0u"W/m^2")
signal = OpticalSignal(lightspot_mask(0.0, 0.0), ls_photocurrent_1uA)

# Photocurrent should be 1uA.
N = 5000
t = 0.5e-3u"s"
L_s5990 = 4.0e-3u"m"
L_s5991 = 9.0e-3u"m"

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
    measured_position_resolution_s5990 / L_s5990 * measurement_photocurrent
expected_noise_s5990 = sqrt(measured_noise_s5990^2 - equivalent_input_current_noise_s5990^2)
resolution_s5990 = L_s5990 * expected_noise_s5990 / measurement_photocurrent

equivalent_input_current_noise_s5991 =
    equivalent_input_voltage_noise / s5991.photodetector.resistance *
    sqrt(measurement_frequency)
measured_position_resolution_s5991 = 1.5e-6u"m"
measured_noise_s5991 =
    measured_position_resolution_s5991 / L_s5991 * measurement_photocurrent
expected_noise_s5991 = sqrt(measured_noise_s5991^2 - equivalent_input_current_noise_s5991^2)
resolution_s5991 = L_s5991 * expected_noise_s5991 / measurement_photocurrent
println("Expected resolution S5990 $resolution_s5990 S5991 $resolution_s5991")

readouts_s5990 = [(SunSensors.read(s5990, signal, t)) / t for _ in 1:N]
readouts_s5991 = [(SunSensors.read(s5991, signal, t)) / t for _ in 1:N]

# Resolution is approximately the standard deviation.
positions_s5990 = [get_pos(s, L_s5990) for s in readouts_s5990]
mean_s5990 = sum.((first.(positions_s5990), last.(positions_s5990))) ./ N
square_errors_s5990 = [(p .- mean_s5990) .^ 2 for p in positions_s5990]
sd_s5990 = sqrt.([sum(first.(square_errors_s5990)), sum(last.(square_errors_s5990))] ./ N)
println(
    "S5990, mean ($(mean_s5990[1]), $(mean_s5990[2])) std ($(sd_s5990[1]), $(sd_s5990[2]))",
)

positions_s5991 = [get_pos(s, L_s5991) for s in readouts_s5991]
mean_s5991 = sum.((first.(positions_s5991), last.(positions_s5991))) ./ N
square_errors_s5991 = [(p .- mean_s5991) .^ 2 for p in positions_s5991]
sd_s5991 = sqrt.([sum(first.(square_errors_s5991)), sum(last.(square_errors_s5991))] ./ N)
println(
    "S5991, mean ($(mean_s5991[1]), $(mean_s5991[2])) std ($(sd_s5991[1]), $(sd_s5991[2]))",
)

photocurrent_s5990 = sum(sum.(readouts_s5990)) / N
photocurrent_s5991 = sum(sum.(readouts_s5991)) / N
println("Photocurrents used in datasheet measurements S5990 1 uA, S5991 1 uA")
println("Photocurrents measured $photocurrent_s5990 $photocurrent_s5991")

## Position detection error
readouts_s5990 = []
readouts_s5991 = []
steps = 100
ls = RayAngles(0.0u"°", 0.0u"°")
for step_x in 1:steps
    for step_y in 1:steps
        # Within 80% of the sensing surface.
        k_x = (-0.5 + step_x / steps) * 0.8
        k_y = (-0.5 + step_y / steps) * 0.8

        x_s5990 = k_x * L_s5990
        y_s5990 = k_y * L_s5990
        s_s5990 =
            SunSensors.read(
                s5990,
                OpticalSignal(lightspot_mask(ustrip(x_s5990), ustrip(y_s5990)), ls),
                t,
            ) / t
        push!(readouts_s5990, ((x_s5990, y_s5990), get_pos(s_s5990, L_s5990)))

        x_s5991 = k_x * L_s5991
        y_s5991 = k_y * L_s5991
        s_s5991 =
            SunSensors.read(
                s5991,
                OpticalSignal(lightspot_mask(ustrip(x_s5991), ustrip(y_s5991)), ls),
                t,
            ) / t
        push!(readouts_s5991, ((x_s5991, y_s5991), get_pos(s_s5991, L_s5991)))
    end
end

error_displacements_s5990 = []
for r in readouts_s5990
    e = r[2] .- r[1]
    push!(error_displacements_s5990, e)
end

errors_s5990 = [sqrt(e[1]^2 + e[2]^2) for e in error_displacements_s5990]

max_error_s5990 = uconvert(u"μm", maximum(errors_s5990))
typical_error_s5990 = uconvert(u"μm", sum(errors_s5990) / length(errors_s5990))
println(
    "Maximum and typical position detection errors in datasheet for S5990 150 μm and 70 μm"
)
println(
    "Maximum and typical position detection errors measured $max_error_s5990 $typical_error_s5990)",
)
n = Int(sqrt(length(errors_s5990)))
h1 = heatmap(reshape(errors_s5990, n, n))

error_displacements_s5991 = []
for r in readouts_s5991
    e = r[2] .- r[1]
    push!(error_displacements_s5991, e)
end

errors_s5991 = [sqrt(e[1]^2 + e[2]^2) for e in error_displacements_s5991]

max_error_s5991 = uconvert(u"μm", maximum(errors_s5991))
typical_error_s5991 = uconvert(u"μm", sum(errors_s5991) / length(errors_s5991))
println(
    "Maximum and typical position detection errors in datasheet for S5991 250 μm and 150 μm"
)
println(
    "Maximum and typical position detection errors $max_error_s5991 $typical_error_s5991)"
)
n = Int(sqrt(length(errors_s5991)))
h2 = heatmap(reshape(errors_s5991, n, n))

plot(h1, h2)
