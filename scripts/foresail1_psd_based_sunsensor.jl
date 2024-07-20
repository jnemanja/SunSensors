include("hamamatsu_s5990.jl")

cover_mask_foresail1_psd = PositionedMask(
    BasicMask(
        [PinholeFeature(0.4e-3u"m")], [Placement(1, Point(0.0, 0.0), 0.0)], 0.2e-3u"m"
    ),
    (0.0, 0.0, 1.33e-3),
    RotX(0.0),
)

# Reading with transimpedance amplifier converter
bits = 10
volts_per_bit = 1.8u"V" / 2^bits
Rf = 9.0e3u"Ω"
bandwidth = 15000.0e6u"Hz"
boltzmann_constant = 1.380649e-23u"J/K"
effective_exposure_time = uconvert(u"s", 1 / bandwidth)
gain = Rf / effective_exposure_time
read_noise = sqrt(4 * boltzmann_constant * s5990.photodetector.temperature * Rf * bandwidth)
adc = ADConverter(gain, volts_per_bit, bits, ustrip(read_noise))

function read_psd(ls::LightSource)
    out = SunSensors.read(s5990, OpticalSignal(cover_mask_foresail1_psd, ls), effective_exposure_time)
    return adc(uconvert.(u"C", out))
end
