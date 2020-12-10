function [] = payload_range_plotter (Velocity, Range, Specific_fuel_consumption, Fuel_volume, Fuel_density, Payload_mass_design, MTOW_design, Crew_members, Maximum_efficiency, a_new, b_new, hold_value)

    hold_value = lower(hold_value);
    m_human = 85;                                                                          % [kg]
    m_crew = Crew_members * m_human;                                                       % [kg]
    SFC = Specific_fuel_consumption;                                                       % [lb/lbh]
    F_V = Fuel_volume;                                                                     % [m³]
    rho = Fuel_density;                                                                    % [kg/m³]
    fuel_mass = F_V * rho;                                                                 % [kg]
    E_max = Maximum_efficiency;
    m_payload_design = Payload_mass_design;                                                % [kg]

    Payload_mass = MTOW_design * (1 - fuel_mass / MTOW_design - a_new * (MTOW_design ^ b_new)) - m_crew;
    Range_a = -log((1 - fuel_mass / MTOW_design) / (0.97 * 0.985 * 0.985 * 0.995)) * Velocity * E_max / SFC;
    empty_mass = (a_new * MTOW_design ^ b_new) * MTOW_design;
    m_payload_null = 0;
    MTOW_empty_withfuel = fuel_mass + empty_mass;
    Range_b = -log((1 - fuel_mass / MTOW_empty_withfuel) / (0.97 * 0.985 * 0.985 * 0.995)) * Velocity * E_max / SFC;
    Y = [m_payload_design; m_payload_design; Payload_mass; m_payload_null];
    X = [m_payload_null; Range; Range_a; Range_b] * 10e-3;

    plot(X, Y, 'Color', rand(1, 3), 'LineWidth', 1.5)
    xlabel('Range [km]')
    ylabel('Payload [kg]')
    title('Payload - Range')
    grid on
    hold on 

    if strcmp(hold_value, "on") == 1

        hold on

    else

        legend('MTOW = 304685 [kg]','MTOW = 268123 [kg]','MTOW = 328146 [kg]','Location','NorthEast')

    end

end