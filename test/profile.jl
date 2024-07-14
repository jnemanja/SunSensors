# Tests Hamamatsu sensor model for S9132.
# The test measures and compares dark current and voltage outputs as reported in the sensors datasheet:
# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s9132_kmpd1075e.pdf

using Statistics
using SunSensors
using Plots
using Meshes

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
