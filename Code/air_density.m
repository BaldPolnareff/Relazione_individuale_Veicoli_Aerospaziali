function [density] = air_density (Sea_level_density, Sea_level_temperature, Altitude, Thermal_gradient)

    density = Sea_level_density * ((Sea_level_temperature + Thermal_gradient * Altitude) / Sea_level_temperature) ^ 4.2561; % [kg/m³]

end