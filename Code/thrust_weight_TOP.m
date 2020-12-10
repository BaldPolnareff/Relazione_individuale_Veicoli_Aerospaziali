function [Thrust_Weight_Ratio, Takeoff_parameter] = thrust_weight_TOP (Takeoff_lift_coeff)

    Takeoff_parameter = 230 * 0.45 / (0.3048 ^ 2);                                      % [kg/m²]
    Thrust_Weight_Ratio = @(x) (x) ./ (Takeoff_parameter * Takeoff_lift_coeff);

end