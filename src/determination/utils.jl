include("psd.jl")
include("profile.jl")

"""
    get_angles(point, positioned_mask)

Calculates the angles of the incident light on the sun sensor from the light spot position given
a positioned mask.
"""
function get_angles(projection_shift::Point, pmask::PositionedMask)
    p = ustrip.(to(projection_shift))
    return rad2deg.(atan.(p / ustrip(pmask.shift[3])))
end

"""
    measure_accuracy_precision(f, point, area_size, samples)

Measures accuracy and precision of the supplied function f in the area_size rectangular
neineighbourhood of the given point, using specified number of samples.

Function f should accept and return a tuple representing the point to be measured at.
"""
function measure_accuracy_precision(
    f::Function, point::Tuple{Float64,Float64}, area_size::Float64, samples::Number
)
    if samples % 2 != 0
        samples = samples - 1
        println("Samples number must be even number. Reducing to $samples.")
    end

    halfsize = area_size / 2
    stepsize = area_size / (samples - 1)
    x_inputs = [point[1] - halfsize + i * stepsize for i in 0:(samples - 1)]
    y_inputs = [point[2] - halfsize + j * stepsize for j in 0:(samples - 1)]

    outs = [f(x, y) .- (x, y) for x in x_inputs for y in y_inputs]
    outs_m = reshape(outs, samples, samples)
    outs_m_x = first.(outs_m)
    outs_m_y = last.(outs_m)

    midpoint = Int(samples / 2)
    x_mean_1 = mean(outs_m_x[:, 1:midpoint])
    x_mean_2 = mean(outs_m_x[:, (midpoint + 1):end])
    y_mean_1 = mean(outs_m_y[1:(midpoint), :])
    y_mean_2 = mean(outs_m_y[(midpoint + 1):end, :])

    k_x = (x_mean_2 - x_mean_1) / midpoint
    k_y = (y_mean_2 - y_mean_1) / midpoint

    biases = [
        (x_mean_1 + k_x * (i - samples / 4), y_mean_1 + k_y * (j - samples / 4)) for
        i in 0:(samples - 1) for j in 0:(samples - 1)
    ]
    debiased = broadcast.(-, outs_m[:], biases)
    std_x = std(first.(debiased))
    std_y = std(last.(debiased))

    return (x_mean_1 + x_mean_2, y_mean_1 + y_mean_2) ./ 2, sqrt(std_x^2 + std_y^2)
end
