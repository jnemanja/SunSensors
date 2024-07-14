const boltzmann_constant = 1.380649e-23u"J/K"
const elementary_charge = 1.602e-19u"C"

"""
    Photodetector

Used for analyzing incident optical signal.
Supports reading intensity and centroid of the optical signal.
"""
abstract type Photodetector end

"""
    charge(photodetector, optical_signal, exposure_time)

Returns the charge of accumulated by the photodetector.
"""
function charge end

"""
    centroid(photodetector, optical_signal, exposure_time)

Returns the centroid of the optical signal received by the photodetector.
"""
function centroid end

"""
    IdealPhotodetector(geometry, responsivity, dark_power_density)

Ideal photodetector model. Reads the optical signal, parametrized with the responsivity factor.
The responsivity factor models the conversion of the received optical signal power to the output
current. Models background and dark current effects with the dark power density.
"""
struct IdealPhotodetector{Dim,T} <: Photodetector
    # Width and height of the photodetector.
    geometry::Polygon{Dim,T}

    # Responsivity factor converts the input optical signal power to the output current intensity.
    # Responsivity is used to represent conversion over the whole spectre of interest.
    responsivity::typeof(1.0u"A/W")

    # Dark power emulates the contribution of the background and dark currents.
    dark_power_density::typeof(1.0u"W/m^2")
end

"""
    IdealPhotodetector(size, responsivity, dark_power_density)

Constructs a rectangular ideal photodetector given a tuple of its width and height,
responsivity factor and dark power.
"""
function IdealPhotodetector(
    size::NTuple{2,typeof(1.0u"m")},
    responsivity::typeof(1.0u"A/W"),
    dark_power_density::typeof(1.0u"W/m^2"),
)
    w = size[1] / 2
    h = size[2] / 2
    geometry = Quadrangle((-w, -h), (w, -h), (w, h), (-w, h))
    return IdealPhotodetector(geometry, responsivity, dark_power_density)
end

"""
    charge(photodetector, no_signal, exposure_time)

Charge accumulation without optical signal present.
"""
function charge(pd::IdealPhotodetector, _::NoSignal, exposure_time::typeof(1.0u"s"))
    return pd.responsivity * area(pd.geometry) * pd.dark_power_density * exposure_time
end

"""
    area_signal(photodetector, optical_signal)

Area of the optical signal over the photodetector's sensing surface.
"""
function area_signal(pd::Photodetector, os::OpticalSignal)
    signals = [intersection(feature, pd.geometry).geom for feature in os.pattern]
    filter!((!isnothing), signals)

    if length(signals) > 0
        area_signal = sum(area.(signals))
    else
        area_signal = 0.0u"m^2"
    end

    return area_signal
end

"""
    charge(ideal_photodetector, optical_signal, exposure_time)

Charge accumulated by the ideal photodetector from the optical signal.
"""
function charge(pd::IdealPhotodetector, os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    power = (
        area(pd.geometry) * pd.dark_power_density + area_signal(pd, os) * os.power_density
    )
    # Here it is assumed that the power contains full emitted spectra and that the
    # responsivity account for the received spectra.
    return pd.responsivity * power * exposure_time
end

function centroid(_::IdealPhotodetector, _::NoSignal, _::typeof(1.0u"s"))
    return nothing
end

"""
    centroid(ideal_photodetector, optical_signal, exposure_time)

Centroid of the optical signal read by the ideal photodetector.
"""
function centroid(pd::IdealPhotodetector, os::OpticalSignal, _::typeof(1.0u"s"))
    signal = []
    for feature in os.pattern
        s = intersection(feature, pd.geometry).geom

        if !isnothing(s)
            push!(signal, s)
        end
    end
    if length(signal) == 0
        return nothing
    end

    cs = to.(Meshes.centroid.(signal))
    as = area.(signal)
    area_signal = sum(as)
    signal_position = sum(cs .* as) / area_signal

    # Shift the position due to the dark current.
    return Point(
        (
            signal_position * area_signal * os.power_density /
            (area_signal * os.power_density + area(pd.geometry) * pd.dark_power_density)
        )...,
    )
end

"""
    NoisyPhotodetector(geometry, responsivity, dark_power_density, temperature, resistance)

Photodetector with shot and thermal noise modeled. Reads the optical signal, parametrized with
the responsivity factor. The responsivity factor models the conversion of the received optical
signal power to the output current. Models background and dark current effects with the dark power
density. Thermal noise is configurable with the temperature and the resistance of the photodetector.
"""
struct NoisyPhotodetector{Dim,T} <: Photodetector
    # Width and height of the photodetector.
    geometry::Polygon{Dim,T}

    # Responsivity factor converts the input optical signal power to the output current intensity.
    # Responsivity is used to represent conversion over the whole spectre of interest.
    responsivity::typeof(1.0u"A/W")

    # Dark power emulates the contribution of the background and dark currents.
    dark_power_density::typeof(1.0u"W/m^2")

    # Photodetector's absolute temperature, used for noise generation.
    # Thermal noise is directly proportional with the temperature.
    temperature::typeof(1.0u"K")

    # Photodetector's resistance, used for noise generation.
    # Thermal noise is inversely proportional with the resistance.
    resistance::typeof(1.0u"Ω")
end

"""
    NoisyPhotodetector(size, responsivity, dark_power_density, temperature, resistance)

Constructs a rectangular noisy photodetector given a tuple of its width and height, responsivity
factor and dark power, and noise parameters of the photodetector's temperature and resistance.
"""
function NoisyPhotodetector(
    size::NTuple{2,typeof(1.0u"m")},
    responsivity::typeof(1.0u"A/W"),
    dark_power_density::typeof(1.0u"W/m^2"),
    temperature::typeof(1.0u"K"),
    resistance::typeof(1.0u"Ω"),
)
    ideal = IdealPhotodetector(size, responsivity, dark_power_density)
    return NoisyPhotodetector(
        ideal.geometry,
        ideal.responsivity,
        ideal.dark_power_density,
        temperature,
        resistance,
    )
end

"""
    standard_deviation(noisy_photodetector, current, exposure_time)

Calculates the standard deviation of the noise for photodetector charge reading.
"""
function standard_deviation(
    pd::NoisyPhotodetector, current::typeof(1.0u"A"), exposure_time::typeof(1.0u"s")
)
    shot_var = 2 * elementary_charge * current / (2 * exposure_time)
    thermal_var =
        4 * boltzmann_constant * pd.temperature / pd.resistance / (2 * exposure_time)
    return sqrt(shot_var + thermal_var)
end

"""
    charge(noisy_photodetector, optical_signal, exposure_time)

Charge accumulated by the noisy photodetector from the optical signal.
"""
function charge(pd::NoisyPhotodetector, os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    c = charge(
        IdealPhotodetector(pd.geometry, pd.responsivity, pd.dark_power_density),
        os,
        exposure_time,
    )
    current = c / exposure_time

    # The shot noise has Poisson distribution. However it is more convenient to use Normal
    # distribution and model several noise effects into a single standard deviation. Poisson is
    # Normal like at high values.
    sd = standard_deviation(pd, current, exposure_time)
    return rand(Normal(ustrip(current), ustrip(sd))) * u"A" * exposure_time
end

"""
    centroid(noisy_photodetector, optical_signal)

Centroid of the optical signal read by the noisy photodetector.
"""
function centroid(pd::NoisyPhotodetector, os::OpticalSignal, exposure_time::typeof(1.0u"s"))
    ipd = IdealPhotodetector(pd.geometry, pd.responsivity, pd.dark_power_density)
    cntr = centroid(ipd, os, exposure_time)
    if isnothing(cntr)
        return nothing
    end

    c = charge(ipd, os, exposure_time)
    current = c / exposure_time
    # Normalize the standard deviation to the photodetector's area.
    sd = standard_deviation(pd, current, exposure_time) / current

    # Scale the standard deviation with the photodetector's geometry.
    bbox = boundingbox(pd.geometry)
    size = bbox.max - bbox.min
    w_sd = ustrip(size[1] * sd)
    h_sd = ustrip(size[2] * sd)

    coord = to(cntr)
    return Point(rand.([Normal(ustrip(coord[1]), w_sd), Normal(ustrip(coord[2]), h_sd)])...)
end
