clear all
clc
close all


    % Variabili fisse e variabili di statistical trend 

        global g = 9.81;                                               % [m/s²]
        global rho_0 = 1.225;                                          % [kg/m³]
        global T_0 = 288;                                              % [K] 
        global z_altitude = 12000;                                     % [m]
        global h_gradient = -6.5e-3;                                   % [k/m]

        M_TO = [250000 ; 275000 ; 251000 ; 351535 ; 352000 ; 227940 ; 365000 ; 242000 ; 380000 ; 230000 ; 158758 ; 362870 ; 447696 ; 299370 ; 78245 ; 273294 ; 103000]; 
        Me = [135500 ; 142400 ; 137000 ; 167829 ; 181400 ; 119950 ; 177755 ; 109400 ; 174000 ; 120600 ; 86069 ; 182480 ; 220128 ; 160530 ; 41413 ; 128808 ; 58300]; 

    % Plot del trend statistico sulla empty mass fraction, con un confronto tra un trend piú grezzo ed uno migliorato (con una formula del Raymer)
    % I valori restituiti ritornano utili dopo.
    
        [a_new, b_new, y_new] = empty_mass_fraction_plotter (M_TO, Me);
       


    %% Studio preliminare

        % Viene scelta la True Air Speed a 10000 m di quota, in quanto un worse case scenario 
        % rispetto al requisito di Mach 0.85

        % Era stato inizialmente preso come punto di design un valore di payload max pari a 50 tonnellate
        % ed un range medio (in crociera) pari a 11000 km, considerando che, diminuendo in parte
        % il payload, si possono raggiungere range piú lunghi.
        % Per problemi di stabilitá del plot nell'intorno di 15000 km, é stato poi scelto un payload di 40 tonnellate
        % per plottare, rispetto al valore di design.

        velocity = 3600 * 0.85 * sqrt(1.4 * 287 * (-50 + 273.15));      % [m/h]
        Range = 11000 * 10 ^ 3;                                         % [m]
        Range_max = 15000 * 10 ^ 3;                                     % [m]
        SFC_1 = 0.478;                                                  % [lb/lbh]
        MTOW_guess = 275000;                                            % [kg]           (Valore del velivolo target (A350), guess iniziale)
        m_payload_design = 50e+3;                                       % [kg]
        m_payload = 40e+3;                                              % [kg]
        N_crew = 8;
        Efficienza_max = 20;

    %% Diagramma MTOW-Range

        [MTOW] = mtow_range_plotter (velocity, Range, N_crew, SFC_1, Efficienza_max, m_payload, MTOW_guess, a_new, b_new);
        mtow_design = MTOW(1);
        mtow_range = MTOW(2);

    %% Diagramma MTOW-Payload 

        % Si considera il range di design, pari a 11000 km

        mtow_payload_plotter (velocity, Range, N_crew, SFC_1, Efficienza_max, MTOW_guess, a_new, b_new);
        
        % Si considera ora il range massimo, pari a 15000 km

        mtow_payload_plotter (velocity, Range_max, N_crew, SFC_1, Efficienza_max, MTOW_guess, a_new, b_new);

    %% Diagramma Payload-Range

        % Si plotta il payload in funzione del range del velivolo, considerando come
        % fuel volume di riferimento quello dell'A350, pari a 141000 L di Avgas.
        % Vengono messi in evidenza i 4 punti notevoli del grafico e i casi notevoli sono
        % riuniti in un plot unico, per avere una comparazione visiva efficace. 

        % Per rappresentare la famiglia di velivoli, se ne considera uno a piú alto range
        % ma minor payload ed uno a piú alto payload ma minor range.
        
        % Rispetto ai valori di design, si ha:

            % Variazione di payload: ± 5'000                            % [kg]
            % Variazione di MTOW:    {0.88*MTOW, 1.077*MTOW}            % [kg]
            % Variazione di Range:   ± 500                              % [km]

        Fuel_volume = 141;                                              % [m³]
        rho_AVGAS = 810;                                                % [kg/m³]
        Range_inf = Range - 500e3;                                      % [m]
        Range_sup = Range + 500e3;                                      % [m]
        Payload_inf = m_payload_design - 5000;                          % [kg]
        Payload_sup = m_payload_design + 5000;                          % [kg]
        MTOW_inf = 0.88 * mtow_design;                                  % [kg]
        MTOW_sup = 1.077 * mtow_design;                                 % [kg]

        figure()
        hold on
        payload_range_plotter (velocity, Range, SFC_1, Fuel_volume, rho_AVGAS, m_payload_design, mtow_design, N_crew, Efficienza_max, a_new, b_new, 'ON');
        payload_range_plotter (velocity, Range_inf, SFC_1, Fuel_volume, rho_AVGAS, Payload_sup, MTOW_sup, N_crew, Efficienza_max, a_new, b_new, 'ON');
        payload_range_plotter (velocity, Range_sup, SFC_1, Fuel_volume, rho_AVGAS, Payload_inf, MTOW_inf, N_crew, Efficienza_max, a_new, b_new, 'off');

        % N.B. Per una questione di implementazione i tre colori delle curve sono casuali e diversi tra loro, potrebbe venire generata una palette
        % non gradita, in tal caso é sufficiente ripetere il plot fino ad uno schema colori ritenuto buono (o editarlo manualmente)

    %% Matching Chart

        % Si fa riferimento ad una quota di 12 km

        rho_air = air_density (rho_0, T_0, z_altitude, h_gradient);     % [kg/m³]

        % Velocitá di stallo

        % Come suggerito dal Raymer, la velocitá di approccio é presa con un valore simile ad altri velivoli
        % della stessa famiglia, dato che non esiste una normativa specifica a riguardo.
        % É stato quindi scelto un valore compreso tra quelli  dell'A350-900 e dell'A350-1000.
        % La velocitá di stallo é ricavata come da normativa.

        velocity_approach = 145 * 0.514444;                             % [m/s] (inizialmente in nodi)
        stall_velocity = velocity_approach / 1.3;                       % [m/s]

        % Una volta consumato il fuel necessario (per il range di 11000 km), 
        % la massa al landing risulta essere circa 214000 kg, valore in linea
        % con i velivoli di riferimento (circa tra 205000 kg e 236000 kg)

        % Secondo quando riportato dal Raymer, valori tipici per il 
        % coefficiente di portanza massimo di un jet transport sono
        % all'incirca tra 2.2 e 3.2.
        % É stato quindi scelto un valore pari a 2.8

        % Prendendo ancora una volta spunto dal Raymer, i calcoli sono svolti
        % con riferimento alla densitá al livello del mare.

        CL_max = 2.8;
        W_S = wing_load_distribution (rho_0, stall_velocity, CL_max);   % [kg/m²]
        figure()
        hold on
        line("xdata", [W_S, W_S], "ydata", [-1000, 1000], "linewidth", 2)
        xlabel('W/S [kg/m²]')
        ylabel('T/W')
        title('Matching Chart') 
        grid on
        
        % Requisiti di crociera

        % N.B. La W/S riportata nei grafici di riferimento esprime la forza
        % in kilogrammi peso, per cui sono state fatte le dovute correzioni dimensionali, 
        % riportando le forze in newton, in modo da ottenere (correttamente)
        % un rapporto spinta-peso adimensionale.

        Velocity = velocity / 3600;                                      % [m/s]
        CD_0 = 0.015;
        Aspect_Ratio = 9;
        e_parameter = 0.921;

        % Rapporto peso-spinta in condizione di velocitá massima

        [T_W_max_speed, max_velocity, K_coeff] = thrust_weight_max_speed (rho_0, Velocity, CD_0, Aspect_Ratio, e_parameter);
        x_max_speed = linspace(0, 700);

        plot(x_max_speed, T_W_max_speed(x_max_speed), 'LineWidth', 1.5)

        % Take Off Distance

        % É stata considerata una ground roll distance di 888 m (circa 3000 ft).
        % La "over 50 ft" distance considerata é maggiore di 1.7 volte.
        % Tale valore é quello a cui corrisponde un TOP (Take Off Parameter)
        % di circa 230 lbs/ft².
        % Inoltre, come suggerito dal Roskam (pag. 107), il coefficiente di 
        % portanza al take off é pari a quello massimo, diviso per un fattore di 1.21

        CL_max_takeoff = 2.5;
        CL_takeoff = CL_max_takeoff / 1.21;

        % Rapporto peso-spinta in base al Take Off Parameter

