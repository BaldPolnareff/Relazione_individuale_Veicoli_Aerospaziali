function [Ratio, Corde_estreme] = wing_profile_plotter (Wing_Surface, Aspect_Ratio, Sweep_Angles, Wing_type)

    Type = lower(Wing_type);

    if strcmp(Type, 'main') == 1

        b = sqrt(Aspect_Ratio * Wing_Surface);
        b_avg = b / 2; 
        c_root = (Wing_Surface / (b_avg) + b_avg * (tand(max(Sweep_Angles)) - tand(min(Sweep_Angles)))) / 2;
        c_tip = c_root - b_avg * (tand(max(Sweep_Angles)) - tand(min(Sweep_Angles)));
        Corde_estreme = [c_root, c_tip];
        Taper_Ratio = c_tip / c_root;
        Ratio = Taper_Ratio;                                                                                        
        y1 = c_root - b_avg * tand(max(Sweep_Angles));
        y2 = y1 - c_tip;
        figure()
        plot([0 b_avg], [0, y2], 'Linewidth', 2)
        hold on
        plot([0 b_avg], [c_root y1], 'Linewidth', 2)
        plot([0 0], [0 c_root], 'Linewidth', 2)
        plot([b_avg b_avg], [y2 y1], 'Linewidth', 2)
        grid on
        title('Main Wing')
        xlabel('X')
        ylabel('Y')

    elseif strcmp(Type, 'tail')  == 1 

        TR = 0.6;
        b = 2 * sqrt(Wing_Surface * (1 - TR) / ((1 + TR) * (tand(max(Sweep_Angles)) - tand(min(Sweep_Angles)))));
        c_root = (b / 2) * (tand(max(Sweep_Angles)) - tand(min(Sweep_Angles))) / (1 - TR);
        c_tip = TR * c_root;
        Corde_estreme = [c_root, c_tip];
        Ratio = TR;
        b_avg = b / 2;
        y1 = c_root - b_avg * tand(max(Sweep_Angles));
        y2 = y1 - c_tip;
        figure()
        plot([0 b_avg], [0, y2], 'Linewidth', 2)
        hold on
        plot([0 b_avg], [c_root y1], 'Linewidth', 2)
        plot([0 0], [0 c_root], 'Linewidth', 2)
        plot([b_avg b_avg], [y2 y1], 'Linewidth', 2)
        grid on
        title('Tail wing')
        xlabel('X')
        ylabel('Y')

    elseif strcmp(Type, 'vertical') == 1

        TR = 0.65;
        AR = 1.72;
        Ratio = TR;
        Wing_Surface = Wing_Surface / 2;
        b = sqrt(AR * Wing_Surface);
        c_root = 2 * Wing_Surface / ((1 + TR) * b);
        c_tip = TR * c_root;
        Corde_estreme = [c_root, c_tip];
        b_avg = b / 2;
        z1 = c_root - b_avg * tand(max(Sweep_Angles));
        z2 = (z1 - c_tip);
        figure()
        plot([0 b_avg], [0, z2], 'Linewidth', 2)
        hold on
        plot([0 b_avg], [c_root z1], 'Linewidth', 2)
        plot([0 0], [0 c_root], 'Linewidth', 2)
        plot([b_avg b_avg], [z2 z1], 'Linewidth', 2)
        grid on
        title('Vertical wing')
        xlabel('X')
        ylabel('Z')

    else

        error('Wrong wing type input, choose either Main, Tail or Vertical')

    end

end