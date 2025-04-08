using SunSensors
using Plots
import Meshes.atol

# Improve precision. Some area miscalculations noticed at the default tolerance.
atol(::Type{Float64}) =  1.0e-12

pinhole_mask = BasicMask(
    [
        PinholeFeature(0.00001u"m"),
        RectangleFeature(0.0001u"m", 0.0001u"m"),
        PinholeFeature(0.00003u"m"),
        RectangleFeature(0.0002u"m", 0.0002u"m"),
    ],
    [
        Placement(1, Point(-0.0002, -0.0006), 0),
        Placement(2, Point(0.0005, -0.0003), 0),
        Placement(3, Point(-0.0005, 0.0003), 0),
        Placement(4, Point(0.0008, 0.0006), 0),
    ],
    0.000005u"m",
)
p_pinhole_mask = PositionedMask(pinhole_mask, (0, 0, 0.0005), RotZ(0))


os = OpticalSignal(p_pinhole_mask, RayAngles(0.000u"°", 0.0u"°"))

##

p = NoisyPhotodetector(
    ((2.0e-3 / 256)u"m", (2.0e-3 / 256)u"m"),
    0.5u"A/W",
    5.0e2u"W/m^2",
    300.0u"K",
    250.0u"Ω",
)
ccd = CcdSensor((256, 256), p)
adc_ccd = ADConverter(1e11u"V/C", 0.2u"V", 8, 0.0)

out1 = SunSensors.read(ccd, os, 0.008u"s")
heatmap(out1)
out2 = adc_ccd(out1)
heatmap(out2)

##

pp = NoisyPhotodetector(
    (2.0e-3u"m" / 256, 2.0e-3u"m"),
    0.5u"A/W",
    1.0e1u"W/m^2",
    300.0u"K",
    20.0u"Ω",
)
prof = ProfileSensor((256, 256), pp)
adc_prof = ADConverter(2e9u"V/C", 0.001u"V", 10, 0.0)

out1 = SunSensors.read(prof, os, 0.0002u"s")
plot(out1[:, 1])
plot!(out1[:, 2])
out2 = adc_prof(out1)
plot(out2[:, 1])
plot!(out2[:, 2])

##

psd = PsdSensor(
    NoisyPhotodetector(
        (2.0e-3u"m", 2.0e-3u"m"),
        0.6u"A/W",
        1.0e-11u"W/m^2",
        300.0u"K",
        50.0u"Ω",
    ),
    [
        2.0e-3; 2.0e-3; 2.0e-3; 2.0e-3;;
        -1; 1; -1; 1;;
        1; 1; -1; -1;;
        -0.0; 0.0; -0.0; 0.0;;
        0.0; -0.0; -0.0; 0.0
    ],
)
adc_psd = ADConverter(1e10u"V/C", 0.0005u"V", 16, 0.1)

out1 = SunSensors.read(psd, os, 0.00002u"s")
println(out1)
out2 = adc_psd(uconvert.(u"C", out1))
println(out2)
