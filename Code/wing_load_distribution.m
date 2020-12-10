function [load] = wing_load_distribution (Density, Stall_velocity, Max_Lift_Coeff)

    g = 9.81;                                                           % [m/s²]
    load = 0.5 * Density * (Stall_velocity ^ 2) * Max_Lift_Coeff / g;   % [kg/m²]

end