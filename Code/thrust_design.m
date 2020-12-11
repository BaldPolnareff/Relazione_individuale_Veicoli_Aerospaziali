function [Thrust] = thrust_design (Thrust_Weight_ratio, MTOW_design, Wing_load_distribution)

    g = 9.81;                                           % [m/s²]
    Thrust = Thrust_Weight_ratio(Wing_load_distribution) * MTOW_design * g;

end