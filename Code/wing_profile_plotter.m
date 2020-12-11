function [Taper_Ratio, Corde_estreme] = wing_profile_plotter (Wing_Surface, Aspect_Ratio, Wing_Surface_Angles, Gamma)

    b = sqrt(Aspect_Ratio * Wing_Surface);
    b_avg = b / 2; 
    c_root = (Wing_Surface / (b_avg) + b_avg * (tand(max(Wing_Surface_Angles)) - tand(min(Wing_Surface_Angles)))) / 2;
    c_tip = c_root - b_avg * (tand(max(Wing_Surface_Angles)) - tand(min(Wing_Surface_Angles)));
    Corde_estreme = [c_root, c_tip];
    Taper_Ratio = c_tip / c_root;
    y1 = c_root - b_avg * tand(max(Wing_Surface_Angles));
    y2 = y1 - c_tip;
    figure()
    plot([0 b_avg], [0, y2])
    hold on
    plot([0 b_avg], [c_root y1])
    plot([0 0], [0 c_root])
    plot([b_avg b_avg], [y2 y1])
    grid on
    xlabel('X')
    ylabel('Y')

end