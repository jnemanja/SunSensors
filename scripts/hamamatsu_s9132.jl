s9132 = ProfileSensor(
    (256, 256),
    NoisyPhotodetector(
        (7.8e-6u"m", 1.9968e-3u"m"),
        0.008u"A/W",
        1.605e-3u"W/m^2",
        298.15u"K",
        5.0e5u"Ω",
    ),
)

low_gain_capacitance = 1e-12u"F"
high_gain_capacitance = 0.2e-12u"F"
adc_low_gain = 1 / low_gain_capacitance
adc_high_gain = 1 / high_gain_capacitance
adc_volt_range = 3.8u"V"
volts_per_bit_8bit = adc_volt_range / 2^8
volts_per_bit_10bit = adc_volt_range / 2^10

adc_8bit_low_gain = ADConverter(adc_low_gain, volts_per_bit_8bit, 8, 0.0)
adc_10bit_low_gain = ADConverter(adc_low_gain, volts_per_bit_10bit, 10, 0.0)
adc_8bit_high_gain = ADConverter(adc_high_gain, volts_per_bit_8bit, 8, 0.0)
adc_10bit_high_gain = ADConverter(adc_high_gain, volts_per_bit_10bit, 10, 0.0)

function s9132_8bit_low_gain(os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    return adc_8bit_low_gain(SunSensors.read(s9132, os, exposure_time))
end

function s9132_10bit_low_gain(os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    return adc_10bit_low_gain(SunSensors.read(s9132, os, exposure_time))
end

function s9132_8bit_high_gain(os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    return adc_8bit_high_gain(SunSensors.read(s9132, os, exposure_time))
end

function s9132_10bit_high_gain(os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    return adc_10bit_high_gain(SunSensors.read(s9132, os, exposure_time))
end
