function [Gross_lift_coeff] = max_gross_lift_coeff_calculator (Mass_coefficients, MTOW_design, Density, Velocity, Stall_velocity, Wing_Surface)

    rho_0 = 1.225;
    g = 9.81;
    coeff = Mass_coefficients;
    k_w = 0.95;
    k_a = 0.4;
    m_cruise_start = coeff(1) * coeff(2) * MTOW_design;
    m_cruise_end = coeff(3) * m_cruise_start;
    average_mass = mean ([m_cruise_start, m_cruise_end]);
    CL_cruise = 2 * average_mass * g / (Density * Velocity ^ 2 * Wing_Surface);
    CL_cruise_wing = CL_cruise / k_w;
    CL_ideal = CL_cruise_wing /k_a;
    CL_max = 2 * MTOW_design * g / (rho_0 * (Stall_velocity * 1.1) ^ 2 * Wing_Surface);
    CL_max_wing = CL_max / k_w;

    Gross_lift_coeff = CL_max_wing / k_a;

end