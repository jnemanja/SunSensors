# Hamamatsu S3931, One-dimensional PSD with 6 mm long sensing surface.

s3931 = PsdSensor(
    NoisyPhotodetector(
        (6.0e-3u"m", 1.0e-3u"m"),
        0.32u"A/W",
        7.85e-5u"W/m^2",
        298.15u"K",
        47.0e3u"Ω",
    ),
    [
        3.0e-3; 3.0e-3;;
        -1.0255; 1.0251;;
        0.0; 0.0;;
        4.05; -4.0;;
        0.0; 0.0
    ],
)
