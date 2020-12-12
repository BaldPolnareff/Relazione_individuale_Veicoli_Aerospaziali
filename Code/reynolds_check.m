function [Re] = reynolds_check (Mean_Aero_Center, Density, Velocity, Dynamic_Viscosity)

    Re = Density * Velocity * Mean_Aero_Center / Dynamic_Viscosity;

    if Re > 500000

        display('Reynolds check: pass');
    
    else 

        error('Reynolds check: fail');

    end


end