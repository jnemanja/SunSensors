# Hamamatsu S5991, Two-dimensional PSD with 9x9 mm sensing surface.

s5991 = PsdSensor(
    NoisyPhotodetector(
        (9.0e-3u"m", 9.0e-3u"m"),
        0.45u"A/W",
        2.96e-5u"W/m^2",
        298.15u"K",
        6.31e3u"Ω",
    ),
    [
        4.5e-3; 4.5e-3; 4.5e-3; 4.5e-3;;
        -1.007; 0.97; 1.04; -1.015;;
        1.008; 1.007; -1.005; -0.982;;
        -7.1; 5.2; -5.8; 8.08;;
        7.5; 8.1; 6.4; 6.1
    ],
)
