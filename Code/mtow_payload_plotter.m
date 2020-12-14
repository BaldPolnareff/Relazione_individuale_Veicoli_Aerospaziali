function [] = mtow_payload_plotter (Velocity, Range, Crew_members, Specific_fuel_comnsumption, Maximum_efficiency, MTOW_guess, a_new, b_new)

    v = Velocity;                                                                       % [m/h]
    E_max = Maximum_efficiency;                                                         % [\]            
    SFC = Specific_fuel_comnsumption;
    R = Range;                                                                          % [m]
    coeff = [0.97; 0.985; (exp(-R * SFC / (v * E_max))); 0.985; 0.995];
    COEFF = prod(coeff);
    m_crew = Crew_members * 85;                                                         % [kg]           (Massa media di 85 kg a persona)
    
    payload = linspace(20e+3, 80e+3, 100);                                              % [kg]

    for i = 1 : (length(payload))
        mtow_p = @(x) x - (m_crew + payload(i)) / (1 - 1.06 * (1 - COEFF) - a_new * (x ^ b_new));
        mtow_payload(i) = fzero(mtow_p, MTOW_guess);
    end

    figure()
    plot(payload, mtow_payload, 'b', 'Linewidth', 2)
    xlabel('payload [kg]')
    ylabel('MTOW [kg]')
    grid on
    
end