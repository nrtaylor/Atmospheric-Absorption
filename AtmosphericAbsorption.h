// NRT: Atmospheric Absorption of Sound

#pragma once

namespace AtmosphericAbsorption
{
    const double kPressureSeaLevelPascals = 101325.0;

    // returns the absorption coefficient in db/meters
    double AbsorptionCoefficient(
        const double frequency_hz,
        const double humidity_percent,
        const double temperature_farenheit,
        const double pressure_pascals = kPressureSeaLevelPascals);

    // returns the gain (0 to 1) 
    double AbsorptionGain(
        const double distance_meters,
        const double frequency_hz,
        const double humidity_percent,
        const double temperature_farenheit,
        const double pressure_pascals = kPressureSeaLevelPascals);

    // returns the frequency in hertz corresponding to the absorption coefficient in db/meters
    // Note: this uses the trignometric solution for cubic equations and assumes:
    // 1. the roots are all real
    // 2. the largest root is the correct soltuion to find the frequency
    // These assumptions have not been proven but have been tested for expected inputs.
    double Frequency(
        const double absorption_coefficient,
        const double humidity_percent,
        const double temperature_farenheit,
        const double pressure_pascals = kPressureSeaLevelPascals);
}