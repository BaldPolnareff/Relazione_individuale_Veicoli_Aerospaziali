function [Passenger_Mass] = passenger_mass_calculator (N_business, N_first_class, N_economy)

    Recommended_mass_each_ec = 82 + 23 + 14;                                      % [kg]
    Recommended_mass_each_other = 82 + 32 + 14;                                 

    Business_mass = N_business * Recommended_mass_each_other;                     % [kg]
    First_class_mass = N_first_class * Recommended_mass_each_other;               % [kg]
    Economy_mass = N_economy * Recommended_mass_each_ec;                          % [kg]

    Passenger_Mass = [Business_mass, First_class_mass, Economy_mass];             % [kg]

end