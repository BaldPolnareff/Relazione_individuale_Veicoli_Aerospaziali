function [Ground_Roll_D, Takeoff_D, Takeoff_D_run, Acc_Go_D, Acc_Stop_D] = meaningful_distances_plotter (N_engines, Thrust, Density, Min_velocity, MTOW_design, friction, Clearance_Height, Aspect_Ratio, C_D_0, Wing_surface)

    g = 9.81;
    Weight = MTOW_design * g;
    e_oswald = 0.92;
    K = 1 / (pi * e_oswald * Aspect_Ratio);
    T = N_engines * Thrust;
    t_liftoff = 2;
    Acceleration = g * (T / Weight - friction);
    Acceleration_oei = g * (Thrust / Weight - friction);
    %b_factor = ((C_D_0 - (friction ^ 2) / 4 * K) * Density / 2 * Weight / Wing_surface);
    b_factor = (Density * g) / (2 * Weight / Wing_surface) * (C_D_0 - (friction ^ 2) / (4 * K));
    V_2 = 1.2 * Min_velocity;
    n_z = (0.734 * V_2 ^ 2) / (Min_velocity ^ 2);
    n_x = 0.3;
    Radius = (V_2 ^ 2) / (g * (n_z - 1));
    gamma_0 = acos(1 - Clearance_Height / Radius);
    d_M = Radius * sin(gamma_0);
    d_liftoff = t_liftoff * V_2;
    Ground_Roll_D = 1 / 2 * b_factor * (log(Acceleration) - log(Acceleration - b_factor * V_2 ^ 2));
    Takeoff_D = Ground_Roll_D + d_liftoff + d_M;
    V_failure = linspace(0, V_2, 10000);
    d_aeo = 1 / (2 * b_factor) * (log(Acceleration) - log(Acceleration - b_factor * (V_failure .^ 2)));
    d_oei= 1 / (2 * b_factor) * (log(Acceleration_oei - b_factor * (V_failure .^ 2)) - log(Acceleration_oei - b_factor * (V_2 ^ 2)));
    Acc_Go_D = d_aeo + d_oei + d_M + d_liftoff;
    d_react = 2 * V_failure;
    d_break = (V_failure .^ 2) / (2 * g * n_x);
    Acc_Stop_D = d_aeo + d_react + d_break;
    Takeoff_D_run = (Ground_Roll_D + d_liftoff + d_M / 2);

    figure()
    plot(V_failure, Acc_Go_D, 'r', 'linewidth', 1.5)
    hold on
    plot(V_failure, Acc_Stop_D, 'b', 'linewidth', 1.5)
    grid on
    title('BFL (Balance Field Length)')
    xlabel('Critical Engine Failure Speed [m/s]')
    ylabel('Distance [m]')
    legend('Acceleration & GO','Acceleration & STOP')

    xc = Ground_Roll_D + d_liftoff;
    yc = Radius;
    xx = linspace(xc, Takeoff_D);
    yy = yc - sqrt(Radius ^ 2 - (xx - xc) .^ 2);
    xx1 = linspace(0, xc);
    yy1 = zeros(1, length(xx1));

    %figure()
    %plot([0 Takeoff_D_run], [0 0], 'LineWidth', 1.5)
    %grid on
    %xlabel('Distance [m]')
    %ylabel('Altitude [m]')
    %title('Take-Off Run')


    %figure()
    %plot([0 Takeoff_D], [0 0], 'LineWidth', 1.5)
    %grid on
    %xlabel('Distance [m]')
    %ylabel('Altitude [m]')
    %title('Take-Off Distance')

    %figure()
    %plot([xx1, xx], [yy1, yy], 'LineWidth', 1.5)
    %grid on
    %xlabel('Distance [m]')
    %ylabel('Altitude [m]')
    %title('Trajectory')

    figure()
    hold on
    plot([xx1, xx], [yy1, yy], 'LineWidth', 3.5)
    plot([0 Takeoff_D], [0 0], 'LineWidth', 3.5)
    plot([0 Takeoff_D_run], [0 0], 'g', 'LineWidth', 3.5)
    grid on
    xlabel('Distance [m]')
    ylabel('Altitude [m]')
    title('Take-Off Distance')
    legend('Trajectory', 'Take-Off Distance', 'Take-Off Run', 'Location', 'NorthWest')

end