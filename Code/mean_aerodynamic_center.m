function [MAC] = mean_aerodynamic_center (Corde_estreme, Taper_Ratio)

    MAC = 2 / 3 * max(Corde_estreme) * (1 + Taper_Ratio + Taper_Ratio ^ 2) / (1 + Taper_Ratio);

end