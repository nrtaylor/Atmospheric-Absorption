// NRT: Atmospheric Absorption of Sound

#pragma once

namespace AtmosphericAbsorption2
{
    const double kPressureSeaLevelPascals = 101325.0;

    // returns the absorption coefficient in db/meters
    double AbsorptionCoefficient(
        const double frequency_hz,
        const double humidity_percent,
        const double temperature_farenheit,
        const double pressure_pascals = kPressureSeaLevelPascals);

    class FilterCutoffSolver {
    public:
        FilterCutoffSolver(const double humidity_percent,
            const double temperature_farenheit,
            const double pressure_pascals = kPressureSeaLevelPascals);

        double Solve(const double distance, const double cutoff_gain = 3.0) const;
    private:
        double nitrogen_relax_freq;
        double oxygen_relax_freq;

        // Precomputed coefficients independent of the absorption coefficient.
        double a1, a2, a3;
    };
}