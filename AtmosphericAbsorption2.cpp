#include "AtmosphericAbsorption2.h"
#include <math.h>
#include <assert.h>

namespace AtmosphericAbsorption2
{
    const double kReferenceAirTemperature = 293.15;    

    double FarenheitToKelvin(const double farenheit)
    {
        return (farenheit + 459.60) * 5 / 9.0;
    }

    // Convert humidity to molar concentration of water vapor as a percentage
    static double HumidityConcentration(
        const double humidity_percent, // 0 to 100.0
        const double temperature_kelvin,
        const double pressure_normalized)
    {
        const double triple_point_temperature_water = 273.16;
        const double csat = -6.8346 * pow(triple_point_temperature_water / temperature_kelvin, 1.261) + 4.6151; // exponent to compute molar concentration
        const double psat = pow(10, csat); // saturation vapor pressure
        return humidity_percent * psat / pressure_normalized; // humidity to molar concentration of water vapor
    }

    static double NitrogenRelaxationFrequency(
        const double humidity_concentration,
        const double temp_normalized,
        double pressure_normalized)
    {
        const double nitrogen_relax_factor = 9 + 280 * humidity_concentration * exp(-4.170 * (pow(temp_normalized, -1.0 / 3) - 1.0));
        return pressure_normalized * (1.0 / sqrt(temp_normalized)) * nitrogen_relax_factor; // an approximate test value is 200
    }

    static double OxygenRelaxationFrequency(
        double humidity_concentration,
        double pressure_normalized)
    {
        const double oxygen_relax_factor = 24 + 40400 * humidity_concentration * (0.02 + humidity_concentration) / (0.391 + humidity_concentration);
        return pressure_normalized * oxygen_relax_factor; // an approximate test value is 25,000
    }

    double AbsorptionCoefficient(
        const double frequency_hz,
        const double humidity_percent,
        const double temperature_farenheit,
        const double pressure_pascals)
    {        
        const double temperature_kelvin = FarenheitToKelvin(temperature_farenheit);        
        const double temp_normalized = temperature_kelvin / kReferenceAirTemperature;
        const double pressure_normalized = pressure_pascals / kPressureSeaLevelPascals;
                                         
        const double humidity_concentration = HumidityConcentration(humidity_percent, temperature_kelvin, pressure_normalized); 

        // Low frequencies are affected more by nitrogen relaxation
        const double nitrogen_relax_freq = NitrogenRelaxationFrequency(humidity_concentration, temp_normalized, pressure_normalized);
        double nitrogen_quantity = 0.1068 * exp(-3352.0 / temperature_kelvin);
        nitrogen_quantity /= (nitrogen_relax_freq + frequency_hz * (frequency_hz / nitrogen_relax_freq));

        // Very high frequencies are affected more by oxygen relaxtion
        const double oxygen_relax_freq = OxygenRelaxationFrequency(humidity_concentration, pressure_normalized);
        double oxygen_quantity = 0.01275 * exp(-2239.10 / temperature_kelvin);
        oxygen_quantity /= (oxygen_relax_freq + (frequency_hz * frequency_hz) / oxygen_relax_freq);

        double pressure_quantity = 1.84e-11 / pressure_normalized;        
        double relaxation_quantity = (nitrogen_quantity + oxygen_quantity) / (temp_normalized * temp_normalized * temp_normalized);

        double absorption_coefficient = 8.686 * sqrt(temp_normalized) * (pressure_quantity + relaxation_quantity);
        absorption_coefficient *= frequency_hz * frequency_hz;        

        return absorption_coefficient;
    }
   
    FilterCutoffSolver::FilterCutoffSolver(const double humidity_percent,
        const double temperature_farenheit, const double pressure_pascals) {

        const double temperature_kelvin = FarenheitToKelvin(temperature_farenheit);
        const double temp_normalized = temperature_kelvin / kReferenceAirTemperature;
        const double pressure_normalized = pressure_pascals / kPressureSeaLevelPascals;

        const double humidity_concentration = HumidityConcentration(humidity_percent,
            temperature_kelvin, pressure_normalized);

        // Low frequencies are affected more by nitrogen relaxation
        nitrogen_relax_freq = NitrogenRelaxationFrequency(
            humidity_concentration, temp_normalized, pressure_normalized);
        // Very high frequencies are affected more by oxygen relaxtion
        oxygen_relax_freq = OxygenRelaxationFrequency(
            humidity_concentration, pressure_normalized);

        const double temp_norm_inv_cube =
            1.0 / (temp_normalized * temp_normalized * temp_normalized);
        const double nitrogen_relax_coefficient =
            temp_norm_inv_cube * 0.1068 * exp(-3352.0 / temperature_kelvin);
        const double oxygen_relax_coefficient =
            temp_norm_inv_cube * 0.01275 * exp(-2239.10 / temperature_kelvin);

        const double pressure_coefficient = (double)1.84e-11 / pressure_normalized;
        // factor multiplied to the absorption quantities
        const double outer_coefficient = 8.686 * sqrt(temp_normalized);

        // Re-arrange the equation as a cubic polynomial with the absorption_coefficient
        // as the constant factor -a4
        // 0 = a1*f^2 + a2*n*f^2/(n^2+f^2) + a3*o*f^2/(o^2+f^2) + a4 where
        // f is the variable frequency, n and o is nitrogen/oxygen relaxation frequencies
        a1 = outer_coefficient * pressure_coefficient;
        a2 = outer_coefficient * nitrogen_relax_coefficient;
        a3 = outer_coefficient * oxygen_relax_coefficient;
    }

    double FindFirstRoot(double a, double b, double c, double d) {
        // Trigonmetric Cubic Solver
        // a, b, c, and d are all real so at least one real root must exist
        // a is > 0.
        // Note: assumes the first root is the largest and correct solution.
        // Convert to depressed cubic using change of variable.
        double p = (3 * a*c - b*b) / (3 * a*a);
        double q = ((2.0 * b*b*b) - (9 * b*a*c) + (27 * a*d*a)) / (27 * a*a*a);

        const double theta = (3 * q * sqrt(-3 / p)) / (2 * p);
        double t0 = 2 * sqrt(-p / 3) * cos(acos(theta) / 3);

        // descriminate = -(4p^3 + 27q^2)

        // For debugging here are the other real roots:
        //const double t2 = -2 * sqrt(-p / 3) * cos(acos(-theta) / 3);
        //const double t1 = -t0 - t2;

        const double root = t0 - b / (3 * a);
        return root;
    }

    double FindRootNewton(double a, double b, double c, double d, double epsilon,
        double x) {
        
        const double discriminant = 18.0*a*b*c*d - 4.0*b*b*b*d + b*b*c*c - 4.0*a*c*c*c - 27.0*a*a*d*d;
        assert(discriminant > 0.00001);

        double a_prime = a * 3;
        double b_prime = b * 2;
        double c_prime = c;

        double x1 = x;
        double delta = epsilon;
        while (fabs(delta) >= epsilon) {
            double cubic = ((a * x1 + b) *  x1 + c) * x1 + d;
            double quadratic = (a_prime * x1 + b_prime) * x1 + c_prime;
            delta = cubic / quadratic;
            x1 = x1 - delta;
        }
        return x1;
    }

    double FilterCutoffSolver::Solve(const double distance, const double cutoff_gain) const {

        const double absorption_coefficient = cutoff_gain / distance;
        const double a4 = -absorption_coefficient;

        const double nitrogen_sq = nitrogen_relax_freq * nitrogen_relax_freq;
        const double oxygen_sq = oxygen_relax_freq * oxygen_relax_freq;

        // expand the denominators (which can then be ignored) and collecting terms
        const double a = a1;
        const double b = a1 * (nitrogen_sq + oxygen_sq) + a2 * nitrogen_relax_freq +
            a3 * oxygen_relax_freq + a4;
        const double c = a1 * nitrogen_sq * oxygen_sq + a2 * oxygen_sq * nitrogen_relax_freq +
            a3 * nitrogen_sq * oxygen_relax_freq + a4 * (nitrogen_sq + oxygen_sq);
        const double d = a4 * oxygen_sq * nitrogen_sq;

        const double root = FindFirstRoot(a, b, c, d);
        const double frequency_hz = sqrt(root);
        return frequency_hz;
    }
}