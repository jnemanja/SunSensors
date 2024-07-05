# Hamamatsu S5990, Two-dimensional PSD with 4x4 mm sensing surface.

s5990 = PsdSensor(
    NoisyPhotodetector(
        (4.0e-3u"m", 4.0e-3u"m"),
        0.45u"A/W",
        6.9e-5u"W/m^2",
        298.15u"K",
        6.055e3u"Ω",
    ),
    [
        2.0e-3; 2.0e-3; 2.0e-3; 2.0e-3;;
        -1.08; 1.02; 1.06; -1.05;;
        1.06; 1.03; -1.009; -1.06;;
        4.1; -2.5; -3.1; 3.3;;
        -2.4; 4.2; -6.1; -4.0
    ],
)
