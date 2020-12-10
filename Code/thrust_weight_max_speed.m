function [Thrust_Weight_Ratio, max_velocity, K_coeff] = thrust_weight_max_speed (Density, Velocity, Drag_Coeff_0, Aspect_Ratio, e_parameter)

    rho_0 = 1.225;                                                                                                                                          % [kg/m³]
    g = 9.81;                                                                                                                                               % [m/s²]
    max_velocity = 1.25 * Velocity;                                                                                                                         % [m/s]
    CD_0 = Drag_Coeff_0;
    K_coeff = 1 / (pi * e_parameter * Aspect_Ratio);
    sigma = Density / rho_0;
    Thrust_Weight_Ratio = @(x) rho_0 * CD_0 * (max_velocity ^ 2) ./ (2 * g * x) + 2 * g * K_coeff * x ./ (Density * sigma * (max_velocity ^ 2));

end