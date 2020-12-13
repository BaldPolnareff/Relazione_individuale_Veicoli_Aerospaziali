function [Avg_Mass] = average_cruise_mass_calculator (MTOW_design, Mass_coefficients)

    % [Beginning, End]
    Mass_extremes = [Mass_coefficients(1) * Mass_coefficients(2) * MTOW_design, Mass_coefficients(1) * Mass_coefficients(2) * Mass_coefficients(3) * MTOW_design];
    Avg_Mass = mean(Mass_extremes);

end