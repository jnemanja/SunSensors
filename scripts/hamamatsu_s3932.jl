# Hamamatsu S3931, One-dimensional PSD with 6 mm long sensing surface.

s3932 = PsdSensor(
    NoisyPhotodetector(
        (12.0e-3u"m", 1.0e-3u"m"),
        0.32u"A/W",
        5.25e-5u"W/m^2",
        298.15u"K",
        89.0e3u"Ω",
    ),
    [
        6.0e-3; 6.0e-3;;
        -1.013; 1.015;;
        0.0; 0.0;;
        1.1; -2.0;;
        0.0; 0.0
    ],
)
