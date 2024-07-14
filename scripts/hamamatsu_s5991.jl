# Hamamatsu S5991, Two-dimensional PSD with 9x9 mm sensing surface.

s5991 = PsdSensor(
    NoisyPhotodetector(
        (9.0e-3u"m", 9.0e-3u"m"),
        0.45u"A/W",
        2.96e-5u"W/m^2",
        298.15u"K",
        7.02e3u"Ω",
    ),
    [
        4.5e-3; 4.5e-3; 4.5e-3; 4.5e-3;;
        -1.006; 0.97; 1.035; -1.015;;
        1.007; 1.007; -1.005; -0.982;;
        -7.1; 5.2; -5.8; 7.58;;
        7.3; 7.8; 6.4; 6.1
    ],
)
