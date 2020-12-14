function [MTOW] = mtow_range_plotter (Velocity, Range, Crew_members, Specific_fuel_comnsumption, Maximum_efficiency, Payload_Mass, MTOW_guess, a_new, b_new)

    v = Velocity;                                                                       % [m/h]
    E_max = Maximum_efficiency;                                                         % [\]            
    SFC = Specific_fuel_comnsumption;
    R = Range;                                                                          % [m]
    coeff = [0.97; 0.985; (exp(-R * SFC / (v * E_max))); 0.985; 0.995];
    COEFF = prod(coeff);
    m_crew = Crew_members * 85;                                                         % [kg]           (Massa media di 85 kg a persona)
    m_payload = Payload_Mass;                                                           % [kg]
    m_to_guess = MTOW_guess;
    m_payload_design = 50e+3;                                                           % [kg]
    range = linspace(10000e+3, 18000e+3, 100);

    m_to = @(x) x - (m_crew + m_payload_design) / (1 - 1.06 * (1 - COEFF) - a_new * (x ^ b_new));
    mtow_design = fzero(m_to, m_to_guess);
    
    for i = 1 : (length(range))
        coeff = [0.97; 0.985; (exp (-range (i) * SFC / (v * E_max))); 0.985; 0.995];
        COEFF = prod(coeff);
        m_to_r = @(x) x - (m_crew + m_payload) / (1 - 1.06 * (1 - COEFF) - a_new * (x ^ b_new));
        mtow_range(i) = fzero(m_to_r, m_to_guess);
    end

    MTOW = [mtow_design, mtow_range];
    figure()
    plot(range, mtow_range, 'Linewidth', 2)
    xlabel('range [m]')
    ylabel('MTOW [kg]')
    grid on

end