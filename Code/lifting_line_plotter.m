function [CL_wing] = lifting_line_plotter (Lift_coeff, Corde_estreme, Wing_span, Aspect_Ratio, Set_Angle, Taper_Ratio)

    i_w = Set_Angle;
    N_segmenti = 200;                                                                                       % [number of segments - 1]
    alpha_twist = -2;                                                                                       % [deg] (Angolo di twist)
    a_2d = 6.3;                                                                                             % [1 / rad] (lift curve slope)
    alpha_0 = -1.5;                                                                                         % [deg] (flap down zero-lift angle of attack)
    theta = pi / (2 * N_segmenti) : pi / (2 * N_segmenti) : pi/2;
    alpha = i_w + alpha_twist : -alpha_twist / (N_segmenti - 1) : i_w;
    
    % Angolo d'attacco del segmento

    z = (Wing_span / 2) * cos(theta);
    MAC_segmento = max(Corde_estreme) * (1 - (1 - Taper_Ratio) * cos(theta));                               % MAC ad ogni segmento
    mu = MAC_segmento * a_2d / (4 * Wing_span);
    LHS = mu .* (alpha - alpha_0) / 57.3;                                                                   % Left Hand Side

    % Risoluzione di N equazioni per trovare i coefficienti A(i)

    for i = 1 : N_segmenti
        for j = 1 : N_segmenti
            B(i, j) = sin((2 * j - 1) * theta(i)) * (1 + (mu(i) * (2 * j - 1)) / sin(theta(i)));
        end
    end
    AA = B \ transpose(LHS);
    for i = 1 : N_segmenti
        sum1(i) = 0;
        sum2(i) = 0;
        for j = 1 : N_segmenti
            sum1(i) = sum1(i) + (2 * j - 1) * AA(j) * sin((2 * j - 1) * theta(i));
            sum2(i) = sum2(i) + AA(j) * sin((2 * j - 1) * theta(i));
        end
    end
    CL = 4 * Wing_span * sum2 ./ MAC_segmento;
    CL_1 = [0, CL];
    y_s = [Wing_span / 2, z];
    figure()
    plot(y_s,CL_1,'-o')
    grid on
    xlabel('Semiala [m]')
    ylabel('Cl - Lift Coefficient')
    title('Andamento Cl - Teoria della Linea Portante')
    CL_wing = pi * Aspect_Ratio * AA(1);                                                                      % Coefficiente di portanza di tutta l'ala

    % Il contributo di portanza dato dall'ala é circa l'85% del totale

    if CL_wing >= 0.85 * Lift_coeff

        display('Il contributo alare alla portanza é sufficiente')

    else 

        error('Wing Lift contribution is below 85%, recheck your values')

    end
    

end