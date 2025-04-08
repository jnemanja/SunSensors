using SunSensors
using Plots
using Meshes

include("../scripts/foresail1_psd_based_sunsensor.jl")
include("../scripts/foresail1_profile_based_sunsensor.jl")

output_dir = "data/output/test/"


## PSD based, sweep and stationary reading plots
steps = 100

angles_sweep = []
measured_angles_sweep = []
measured_angles_stationary = []
for i in 0:(steps - 1)
    angle_x = (-0.5 + i / (steps - 1)) * 120.0u"°"
    signal = RayAngles(angle_x, 0.0u"°")
    out_sweep = read_psd(signal)
    a = get_angles(Point(spot_position(out_sweep, s5990)...), cover_mask_foresail1_psd)
    push!(angles_sweep, (angle_x, 0.0u"°"))
    push!(measured_angles_sweep, a)

    out_stationary = read_psd(RayAngles(0.0u"°", 0.0u"°"))
    a = get_angles(Point(spot_position(out_stationary, s5990)...), cover_mask_foresail1_psd)
    push!(measured_angles_stationary, a)
end

plot(
    first.(angles_sweep),
    first.(measured_angles_sweep);
    xlabel="light ray angle",
    ylabel="measured angle (°)",
)
savefig(output_dir * "foresail1_psd_based_sunsensor_x_sweep.svg")

plot(first.(measured_angles_stationary); xlabel="sample", ylabel="measured angle (°)")
savefig(output_dir * "foresail1_psd_based_sunsensor_x_stationary.svg")

println("Stationary and sweep plots of the PSD based sunsensor model saved in $output_dir")

## PSD based, accuracy and precision
steps = 20
samples = 20
function measure(x, y)
    signal = RayAngles(x * u"°", y * u"°")
    out = read_psd(signal)
    angles = get_angles(Point(spot_position(out, s5990)...), cover_mask_foresail1_psd)
    return angles
end

angles = []
accuracies = []
precisions = []
for i in 0:(steps - 1)
    for j in 0:(steps - 1)
        angle_x = (-0.5 + i / (steps - 1)) * 110u"°"
        angle_y = (-0.5 + j / (steps - 1)) * 110u"°"
        a, p = measure_accuracy_precision(
            measure, ustrip.((angle_x, angle_y)), 0.05e-3, samples
        )
        push!(angles, (angle_x, angle_y))
        push!(accuracies, a)
        push!(precisions, p)
    end
end

accuracy_psd = reshape(sqrt.(sum.(broadcast.(x -> x .^ 2, accuracies))), steps, steps)

angles_axes = [(-0.5 + i / (steps - 1)) * 110u"°" for i in 0:(steps - 1)]
heatmap(
    angles_axes,
    angles_axes,
    accuracy_psd;
    label="ASS angle detection error",
    xlabel="x angle",
    ylabel="y angle",
    colorbar_title="accuracy (°)",
)
savefig(output_dir * "foresail1_psd_based_sunsensor_accuracy.svg")

heatmap(
    angles_axes,
    angles_axes,
    reshape(precisions, steps, steps);
    label="ASS angle precision",
    xlabel="x angle",
    ylabel="y angle",
    colorbar_title="precision (°)",
)
savefig(output_dir * "foresail1_psd_based_sunsensor_precision.svg")

println("Accuracy and precision plots of the PSD based sunsensor model saved in $output_dir")

## Profile based, sweep and stationary reading plots
steps = 100

angles_sweep = []
measured_angles_sweep = []
outs_sweep = []
measured_angles_stationary = []
outs_stationary = []
for i in 0:(steps - 1)
    angle_x = (-0.5 + i / (steps - 1)) * 80.0u"°"
    signal = RayAngles(angle_x, 0.0u"°", 150.0u"W/m^2")
    out_sweep = read_profile(signal, 0.001u"s")
    a = get_angles(Point(spot_position(out_sweep, s9132)...), cover_mask_foresail1_profile)
    push!(outs_sweep, out_sweep)
    push!(angles_sweep, (angle_x, 0.0u"°"))
    push!(measured_angles_sweep, a)

    out_stationary = read_profile(RayAngles(0.0u"°", 0.0u"°", 150.0u"W/m^2"), 0.001u"s")
    a = get_angles(
        Point(spot_position(out_stationary, s9132)...), cover_mask_foresail1_profile
    )
    push!(outs_stationary, out_stationary)
    push!(measured_angles_stationary, a)
end

plot(
    first.(angles_sweep),
    first.(measured_angles_sweep);
    xlabel="light ray angle",
    ylabel="measured angle (°)",
)
savefig("data/output/test/foresail1_profile_based_sunsensor_x_sweep.svg")

plot(first.(measured_angles_stationary); xlabel="sample", ylabel="measured angle (°)")
savefig("data/output/test/foresail1_profile_based_sunsensor_x_stationary.svg")

outs_sweep_x = reshape(vcat([o[:, 1] for o in outs_sweep]...), 256, steps)
heatmap(outs_sweep_x; xlabel="sample", ylabel="pixel readings")
savefig("data/output/test/foresail1_profile_based_sunsensor_x_sweep_raw.svg")

outs_stationary_x = reshape(vcat([o[:, 1] for o in outs_stationary]...), 256, steps)
heatmap(outs_stationary_x; xlabel="sample", ylabel="pixel readings")
savefig("data/output/test/foresail1_profile_based_sunsensor_x_stationary_raw.svg")

println("Stationary and sweep plots of the PSD based sunsensor model saved in $output_dir")

## Profile based, accuracy and precision
steps = 20
samples = 20
function measure(x, y)
    signal = RayAngles(x * u"°", y * u"°")
    try
        out = read_profile(signal, 0.0001u"s")
        return get_angles(Point(spot_position(out, s9132)...), cover_mask_foresail1_profile)
    catch _
        println("Exception: $x $y")
        return (x, y)
    end
end

angles = []
accuracies = []
precisions = []
for i in 0:(steps - 1)
    for j in 0:(steps - 1)
        angle_x = (-0.5 + i / (steps - 1)) * 36u"°"
        angle_y = (-0.5 + j / (steps - 1)) * 36u"°"
        a, p = measure_accuracy_precision(
            measure, ustrip.((angle_x, angle_y)), 0.5e-3, samples
        )
        push!(angles, (angle_x, angle_y))
        push!(accuracies, a)
        push!(precisions, p)
    end
end

accuracy_profile = reshape(sqrt.(sum.(broadcast.(x -> x .^ 2, accuracies))), steps, steps)
angles_axes = [(-0.5 + i / (steps - 1)) * 36u"°" for i in 0:(steps - 1)]

heatmap(
    angles_axes,
    angles_axes,
    accuracy_profile;
    label="DSS angle detection error",
    xlabel="x angle",
    ylabel="y angle",
    colorbar_title="accuracy (°)",
)
output_dir = "data/output/test/"
savefig(output_dir * "foresail1_profile_based_sunsensor_accuracy.svg")

heatmap(
    angles_axes,
    angles_axes,
    reshape(precisions, steps, steps);
    label="DSS angle precision",
    xlabel="x angle",
    ylabel="y angle",
    colorbar_title="precision (°)",
)
savefig(output_dir * "foresail1_profile_based_sunsensor_precision.svg")

println("Accuracy and precision plots of the profile based sunsensor model saved in $output_dir")
