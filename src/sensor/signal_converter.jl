"""
    SignalConverter

Models conversions of the analog signal read from the photodetectors.
"""
abstract type SignalConverter end

"""
    convert(signal_converter, analog_in)

Converts the analog input signal using the signal converter model.
"""
function convert end

"""
    ADConverter(gain, volts_per_bit, bits, read_noise)

Models conversion of the analog signal read from the photodetector to digital signal.
The converter is configurable with the signal gain, volts per bit factor, the number of the bits
for the output signal and the read noise standard deviation.
"""
struct ADConverter <: SignalConverter
    # Conversion gain, charge to volts.
    gain::typeof(1.0u"V/C")

    # Volts per bit.
    volts_per_bit::typeof(1.0u"V")

    # Number of bits for the output signal.
    bits::UInt8

    # Read noise standard deviation.
    read_noise::Float64
end

"""
    convert(ad_converter, analog_in)

Converts the analog input signal to digital using the analog to digital converter model.
"""
function convert(adc::ADConverter, analog_in::typeof(1.0u"C"))
    signal = rand(Normal(ustrip(analog_in * adc.gain), adc.read_noise)) * 1.0u"V"
    return clamp(round(signal / adc.volts_per_bit), 0, 2^adc.bits - 1)
end

"""
    AAConverter(gain, read_noise)

Models conversion of the analog signal read from the photodetector to analog signal.
The converter is configurable with the signal gain and the read noise standard deviation.
"""
struct AAConverter <: SignalConverter
    # Conversion gain, charge to volts.
    gain::typeof(1.0u"V/C")

    # Read noise standard deviation.
    read_noise::Float64
end

"""
    convert(aa_converter, analog_in)

Converts the analog input signal to analog using the analog to analog converter model.
"""
function convert(aac::AAConverter, analog_in::typeof(1.0u"C"))
    return rand(Normal(analog_in * aac.gain, aac.read_noise))
end

(sc::SignalConverter)(input_signal::typeof(1.0u"C")) = convert(sc, input_signal)
function (sc::SignalConverter)(input_signal::Vector{typeof(1.0u"C")})
    return map(s -> convert(sc, s), input_signal)
end
function (sc::SignalConverter)(input_signal::Matrix{typeof(1.0u"C")})
    return map(s -> convert(sc, s), input_signal)
end
