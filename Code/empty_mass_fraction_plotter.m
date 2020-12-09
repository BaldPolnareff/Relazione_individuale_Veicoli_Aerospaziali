function [a_new, b_new, y_new] = empty_mass_fraction_plotter (Takeoff_mass, Empty_mass)

    x = 1e4:1e6;
    A = 0.97;
    b = -0.06;
    y_lett = A .* (x .^ b);
    x_1 = Takeoff_mass;
    y_1 = Empty_mass ./ Takeoff_mass;
    
    figure()
    semilogx(x, y_lett, '-', x_1(1), y_1(1), 'o', x_1(2), y_1(2), 'o', x_1(3), y_1(3), 'o', x_1(4), y_1(4), 'o', x_1(5), y_1(5), 'o', x_1(6), y_1(6), 'o', x_1(7), y_1(7), 'o', x_1(8), y_1(8), 'o', x_1(9), y_1(9), 'o', x_1(10), y_1(10), 'o', x_1(11), y_1(11), 'o', x_1(12), y_1(12), 'o', x_1(13), y_1(13), 'o', x_1(14), y_1(14), 'o', x_1(15), y_1(15), 'o', x_1(16), y_1(16), 'o')
    title("Empty mass fraction trend")
    xlabel('Takeoff Weight [kg]')
    ylabel('Empty weight fraction')
    grid on

    % Trend migliorato con formula del Raymer

    n = length(x_1);
    A1 = [(ones(n, 1)) log(x_1)];
    c = A1 \ log(y_1);
    a_new = exp(c(1));
    b_new = c(2);
    y_new = a_new .* (x .^ b_new);
    hold on
    semilogx(x, y_new, 'r--')
    legend('Statistical trend', 'B-787', 'A350', 'A330 Neo-900', 'B777-300', 'B777 X-900', 'B787-8', 'A340-500', 'A330-3000', 'A340-600', 'A330-200', 'B767-300', 'B747-400', 'B747-8', 'B777-300', 'B737-800', 'MD11', 'Improved trend', 'location', 'SouthWest')

end