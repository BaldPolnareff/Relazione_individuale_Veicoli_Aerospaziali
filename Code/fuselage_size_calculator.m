function [Business_dim, FC_dim, Economy_dim, Fuselage_Total_Length] = fuselage_size_calculator (Business_vec, FC_vec, Economy_vec, N_seats)

    Toilet = [4, 1];                                                                                        % [/, m] (Number, length)

    % BUSINESS 1-2-1
    width_bn = 4 * Business_vec(3) + 2 * Business_vec(1);                                                   % [m]
    length_bn = N_seats(1) / 4 * Business_vec(2);                                                           % [m]

    % FIRST CLASS 2-2-2
    width_fc = 6 * FC_vec(3) + 2 * FC_vec(1);                                                               % [m]
    length_fc = N_seats(2) / 6 * FC_vec(2);                                                           % [m]

    % ECONOMY 3-3-3
    width_ec = 9 * Economy_vec(3) + 2 * Economy_vec(1);                                                     % [m]
    length_ec = N_seats(3) / 9 * Economy_vec(2);                                                            % [m] 

    Business_dim = [width_bn, length_bn];                                                                   % [m, m]
    FC_dim = [width_fc, length_fc];                                                                         % [m, m]
    Economy_dim = [width_ec, length_ec];                                                                    % [m, m]

    Extra_lengths = [(1.75 * (Economy_dim(1) + 2 * 0.102)),(4 * 1.07), (Toilet(1) * Toilet(2)), (1 + 1.75 * (Economy_dim(1) + 2 * 0.102))];
    % [nose, portelloni, bagni, coda]

    Fuselage_Total_Length = sum(Extra_lengths) + Business_dim(2) + FC_dim(2) + Economy_dim(2);

end