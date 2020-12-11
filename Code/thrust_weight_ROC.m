function [Thrust_Weight_Ratio] = thrust_weight_ROC (Rate_of_climb, Drag_coeff_takeoff_0, K_coeff, Max_Efficiency)

    g = 9.81;                                                                                                                           % [m/s²]
    rho_0 = 1.225;                                                                                                                      % [kg/m³]
    Thrust_Weight_Ratio = @(x) Rate_of_climb ./ sqrt(2 * g * x ./ (rho_0 * sqrt(Drag_coeff_takeoff_0 / K_coeff))) + 1 / Max_Efficiency;

end