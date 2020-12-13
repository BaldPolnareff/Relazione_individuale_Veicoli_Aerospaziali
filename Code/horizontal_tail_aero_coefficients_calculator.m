function [C_M0_Wing_Body, Tail_Lift_Coeff] = horizontal_tail_aero_coefficients_calculator (C_MO_af, Incidence, Aspect_Ratio, Sweep_Angles, Cruise_Lift_coeff, Volume_coeff)

    C_M0_Wing_Body = C_MO_af * (Aspect_Ratio * cosd(max(Sweep_Angles)) ^ 2) / (Aspect_Ratio + 2 * cosd(max(Sweep_Angles))) + 0.01 * Incidence;
    Tail_Lift_Coeff = (C_M0_Wing_Body + Cruise_Lift_coeff * 0.2) / Volume_coeff;

end