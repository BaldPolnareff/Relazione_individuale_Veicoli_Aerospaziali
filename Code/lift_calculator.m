function [Lift] = lift_calculator (Density, Velocity, Wing_surface, Lift_coeff)

    Lift = 0.5 * Density * Velocity ^ 2 * Wing_surface * Lift_coeff;

end