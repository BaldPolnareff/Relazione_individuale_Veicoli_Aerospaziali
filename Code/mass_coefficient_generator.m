function [coeff, COEFF] = mass_coefficient_generator (Specific_fuel_consumption, Max_Efficiency, Velocity, Range)

    range = linspace(10000e+3, 18000e+3, 100);
    coeff = [0.97; 0.985; (exp(-Range * Specific_fuel_consumption / (Velocity * Max_Efficiency))); 0.985; 0.995];
    COEFF = prod(coeff);

end