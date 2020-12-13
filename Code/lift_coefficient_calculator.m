function [Lift_coeff] = lift_coefficient_calculator (Average_mass, Density, Velocity, Wing_surface)

    g = 9.81;
    Lift_coeff = 2 * Average_mass * g / (Density * Velocity ^ 2 * Wing_surface);

end