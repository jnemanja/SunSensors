# Tests Hamamatsu sensor model for S9132.
# The test measures and compares dark current and voltage outputs as reported in the sensors datasheet:
# https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s9132_kmpd1075e.pdf

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

global current = 0.0u"C*s^-1"
global bits_8_low = 0
global bits_10_low = 0
global bits_8_high = 0
global bits_10_high = 0
samples = 5000
N = 512
readings = []
for _ in 1:samples
    push!(readings, SunSensors.read(s9132, signal, t))
end

dark_current = sum(sum.(readings)) / (t * N * samples)

bits_8_low = sum(sum.(adc_8bit_low_gain.(readings))) / (N * samples)
bits_10_low = sum(sum.(adc_10bit_low_gain.(readings))) / (N * samples)
bits_8_high = sum(sum.(adc_8bit_high_gain.(readings))) / (N * samples)
bits_10_high = sum(sum.(adc_10bit_high_gain.(readings))) / (N * samples)
println("Dark current in datasheet 0.2 pA")
println("Dark current $dark_current")
println("Dark output voltages in datasheet are 20 mV for low gain and 100 mV for high gain")
println("Converted to output bits, for low gain 8-bit is 1, 10-bit is 5 and for high gain 8-bit is 7 and 10-bit is 27")
println("Dark bits measured $bits_8_low, $bits_10_low, $bits_8_high, $bits_10_high")
println("Noise significantly biases the outputs, given the signal is clammped above 0")
