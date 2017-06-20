#include "AtmosphericAbsorption.h"
#include <math.h>

namespace AtmosphericAbsorption
{
    const double kCelsius20ToKelvin = 293.15;    

    static double FarenheitToKelvin(
        const double farenheit)
    {
        return (farenheit + 459.60) * 5 / 9.0;
    }

    // humidity to molar concentration of water vapor
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
        const double temp_normalized = temperature_kelvin / kCelsius20ToKelvin;
        const double pressure_normalized = pressure_pascals / kPressureSeaLevelPascals;
                                         
        const double humidity_concentration = HumidityConcentration(humidity_percent, temperature_kelvin, pressure_normalized); 

        // Low frequencies are affected more by nitrogen relaxation
        const double nitrogen_relax_freq = NitrogenRelaxationFrequency(humidity_concentration, temp_normalized, pressure_normalized);
        double nitrogen_quantity = 0.1068 * exp(-3352.0 / temperature_kelvin);
        nitrogen_quantity /= (nitrogen_relax_freq + frequency_hz * (frequency_hz / nitrogen_relax_freq));

        // Very high frequencies are affect more by oxygen relaxtion
        const double oxygen_relax_freq = OxygenRelaxationFrequency(humidity_concentration, pressure_normalized);
        double oxygen_quantity = 0.01275 * exp(-2239.10 / temperature_kelvin);
        oxygen_quantity /= (oxygen_relax_freq + (frequency_hz * frequency_hz) / oxygen_relax_freq);

        double pressure_quantity = 1.84e-11 / pressure_normalized;        
        double relaxation_quantity = (nitrogen_quantity + oxygen_quantity) / (temp_normalized * temp_normalized * temp_normalized);

        double absorption_coefficient = 8.686 * sqrt(temp_normalized) * (pressure_quantity + relaxation_quantity);
        absorption_coefficient *= frequency_hz * frequency_hz;        

        return absorption_coefficient;
    }
    
    double AbsorptionGain(
        const double distance_meters,
        const double frequency_hz,
        const double humidity_percent,
        const double temperature_farenheit,
        const double pressure_pascals)
    {
        const double absorption_coefficient = AbsorptionCoefficient(frequency_hz, humidity_percent, temperature_farenheit, pressure_pascals);
        const double absorption_factor = pow(10.0, -(absorption_coefficient * distance_meters) / 20);
        return absorption_factor;
    }

    double Frequency(
        double absorption_coefficient,
        double humidity_percent,
        double temperature_farenheit,
        double pressure_pascals)
    {
        const double temperature_kelvin = FarenheitToKelvin(temperature_farenheit);        
        const double temp_normalized = temperature_kelvin / kCelsius20ToKelvin;
        const double pressure_normalized = pressure_pascals / kPressureSeaLevelPascals;

        const double humidity_concentration = HumidityConcentration(humidity_percent, temperature_kelvin, pressure_normalized);

        const double nitrogen_relax_freq = NitrogenRelaxationFrequency(humidity_concentration, temp_normalized, pressure_normalized);
        const double oxygen_relax_freq = OxygenRelaxationFrequency(humidity_concentration, pressure_normalized);

        const double temp_norm_inv_cube = 1.0 / (temp_normalized * temp_normalized * temp_normalized);
        const double nitrogen_relax_coefficient = temp_norm_inv_cube * 0.1068 * exp(-3352.0 / temperature_kelvin);
        const double oxygen_relax_coefficient = temp_norm_inv_cube * 0.01275 * exp(-2239.10 / temperature_kelvin);

        const double pressure_coefficient = (double)1.84e-11 / pressure_normalized;
        const double outer_coefficient = 8.686 * sqrt(temp_normalized); // factor multiplied to the absorption quantities

        // re-arrange the equation as a cubic polynomial with the absorption_coefficient as the constant factor -a4
        // 0 = a1*f^2 + a2*n*f^2/(n^2+f^2) + a3*o*f^2/(o^2+f^2) + a4 where
        // f is the variable frequency, n and o is nitrogen/oxygen relaxation frequencies
        const double a1 = outer_coefficient * pressure_coefficient;
        const double a2 = outer_coefficient * nitrogen_relax_coefficient;
        const double a3 = outer_coefficient * oxygen_relax_coefficient;
        const double a4 = -absorption_coefficient;

        const double nitrogen_sq = nitrogen_relax_freq * nitrogen_relax_freq;
        const double oxygen_sq = oxygen_relax_freq * oxygen_relax_freq;

        // expand the denominators (which can then be ignored) and collecting terms
        const double a = a1;
        const double b = a1 * nitrogen_sq + a1 * oxygen_sq + a2 * nitrogen_relax_freq + a3 * oxygen_relax_freq + a4;
        const double c = a1 * nitrogen_sq * oxygen_sq + a2 * oxygen_sq * nitrogen_relax_freq + a3 * nitrogen_sq * oxygen_relax_freq + a4 * nitrogen_sq + a4 * oxygen_sq;
        const double d = a4 * oxygen_sq * nitrogen_sq;

        // Trigonmetric Cubic Solver
        // Note: assumes the first root is the largest and correct solution.
        double p = (3 * a*c - b*b) / (3 * a*a);
        double q = ((2.0 * b*b*b) - (9 * b*a*c) + (27 * a*d*a)) / (27 * a*a*a);

        const double theta = (3 * q * sqrt(-3 / p)) / (2 * p);
        double t0 = 2 * sqrt(-p / 3) * cos(acos(theta) / 3);
        
        // For debugging here are the other real roots:
        //const double t2 = -2 * sqrt(-p / 3) * cos(acos(-theta) / 3);
        //const double t1 = -t0 - t2;

        t0 -= b / (3 * a);

        const double frequency_hz = sqrt(t0);
        return frequency_hz;
    }
}