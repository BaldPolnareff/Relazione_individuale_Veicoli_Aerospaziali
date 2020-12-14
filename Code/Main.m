clear all
clc
close all


    % Variabili fisse e variabili di statistical trend 

        g = 9.81;                                               % [m/s²]
        rho_0 = 1.225;                                          % [kg/m³]
        T_0 = 288;                                              % [K] 
        z_altitude = 12000;                                     % [m]
        h_gradient = -6.5e-3;                                   % [k/m]

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
        sigma_0 = 1;
        W_S = wing_load_distribution (rho_0, stall_velocity, CL_max);   % [kg/m²]
        figure()
        hold on
        line("xdata", [W_S, W_S], "ydata", [-1000, 1000], "linewidth", 2)
        hold on 
        line("xdata", [-1000, 1000], "ydata", [1 / (sigma_0 * Efficienza_max), 1 / (sigma_0 * Efficienza_max)], "linewidth", 3)
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
        x_space = linspace(0, 700);

        plot(x_space, T_W_max_speed(x_space), 'LineWidth', 1.5)

        % Take Off Distance

        % É stata considerata una ground roll distance di 888 m (circa 3000 ft).
        % La "over 50 ft" distance considerata é maggiore di 1.7 volte.
        % Tale valore é quello a cui corrisponde un TOP (Take Off Parameter)
        % di circa 230 lbs/ft².
        % Inoltre, come suggerito dal Roskam (pag. 107), il coefficiente di 
        % portanza al take off é pari a quello massimo, diviso per un fattore di 1.21

        CL_max_takeoff = 2.5;
        CL_takeoff = CL_max_takeoff / 1.21;
        CD_0_takeoff = 0.075 + 2;
        ROC = 3000 * 0.00508;                                               % [m/s]

        % Rapporto peso-spinta in base al Take Off Parameter

        [T_W_takeoff, TOP] = thrust_weight_TOP (CL_takeoff);

        plot(x_space, T_W_takeoff(x_space, 'Linewidth', 1.5))
        axis([0 700 0 1])

        % Rapporto peso-spinta in base al Rate of Climb

        T_W_roc = thrust_weight_ROC (ROC, CD_0_takeoff, K_coeff, Efficienza_max);
        plot(x_space, T_W_roc(x_space), 'LineWidth', 1.5)

        % Punto di design

        % Il punto di design scelto é il punto in grado di minimizzare la superficie alare, 
        % a discapito di un leggero aumento della spinta necessaria.
        % Tale scelta é motivata dalla necessitá di ridurre la superficie alare,
        % inizialmente troppo alta. Inoltre, il motore scelto ha una spinta
        % disponibile piú che sufficiente a coprire la spinta necessaria,
        % anche con un ampio margine.

        plot(W_S, T_W_takeoff(W_S), 'or', 'Linewidth', 3)
        legend('STALLO', 'CRUISE', 'TOP', 'ROC', 'CEILING', 'Design Point', 'Location', 'NorthEast')

%% Progetto Ala

        % Spinta di design

        % N.B. Si fa riferimento ad una massa media tra inizio e fine crociera

        % % Secondo indicazione del Sadraey, é stata ipotizzata la presenza di 
        % HLD (flap e slat) con un coefficiente di portanza rispettivamente di
        % 0.9 e 0.4 (pag. 236), dei quali si tiene conto nel calcolare la portanza.

        % Tra i profili segnalati era stato inzialmente scelto il GOE 682, con uno sweep angle di
        % 30 deg, come suggerito dal Raymer (pag. 53).
        % Ci si é poi resi conto che la scelta iniziale non era ottimale, poiché avrebbe
        % fornito un calettamento di circa 1 deg, mentre i valori tipici sono compresi tra 3 e 5 deg.
        % É stato quindi scelto il profilo NACA1412, poiché fornisce un coefficiente di portanza
        % ideale per 3.5 deg, oltre che possededere un coefficiente di portanza massimo di 1.51,
        % a fronte del valore richiesto di 1.3
        % Questo valore é abbastanza conservativo e fornisce un buon margine di sicurezza dallo
        % Il Reynolds massimo riportato da AirfoilTools é di 10^6. 
        % Una volta calcolata la corda media aerodinamica, si potrá
        % verificare il Reynolds ottenuto e validare (o scartare) i dati presi dal sito.
     
        [coeff, COEFF] = mass_coefficient_generator (SFC_1, Efficienza_max, velocity, Range);   % (coefficienti di massa in coeff, loro produttoria in COEFF)

        Thrust_design = thrust_design (T_W_takeoff, mtow_design, W_S);
        Surface_w = mtow_design / W_S;                                                          % [m²]

        CL_max_gross = max_gross_lift_coeff_calculator (coeff, mtow_design, rho_air, Velocity, stall_velocity, Surface_w);

        % Dimensionamento geometrico dell'ala

        % La corda media aerodinamica viene calcolata con una formula approssimata,
        % fornita dal Raymer (pag. 104) ed in seguito usata anche per gli impennaggi

        Sweep_Angles = [29 18];                                                                 % [deg]
        Gamma = 5;                                                                              % [deg]

        [Taper_Ratio, Corde_estreme] = wing_profile_plotter(Surface_w, Aspect_Ratio, Sweep_Angles, 'Main');
        MAC = mean_aerodynamic_center(Corde_estreme, Taper_Ratio);
        Dynamic_Viscosity = 1.46e-5;
        Reynolds = reynolds_check(MAC, rho_air, Velocity, Dynamic_Viscosity);

%% Progetto Fusoliera

        % I parametri sono stati scelti seguendo i suggerimenti del Sadraey (pag. 376)

        N_business = 26;                                                                        % Posti in business
        N_first_class = 48;                                                                     % Posti in first class
        N_economy = 315;                                                                        % Posti in economy (9 per fila)

        Passengers_mass = passenger_mass_calculator(N_business, N_first_class, N_economy);      % [kg kg kg] (Massa Business, Massa First Class, Massa Economy)

        N_seats = [N_business, N_first_class, N_economy];                                       
        Business_geometry = [70, 104, 75] * 10 ^ -2;                                            % [m, m, m] (Larghezza corridoio, Spazio tra sedili, Larghezza sedili)
        First_class_geometry = [60, 92, 60] * 10 ^ -2;                                          % [m, m, m] (Larghezza corridoio, Spazio tra sedili, Larghezza sedili)
        Economy_geometry = [50, 72, 46] * 10 ^ -2;                                              % [m, m, m] (Larghezza corridoio, Spazio tra sedili, Larghezza sedili)

        % Configurazione scelta:

        % Business:     1-2-1
        % First Class:  2-2-2 
        % Economy:      3-3-3

        % Alle dimensioni é stato preventivamente aggiunto uno spessore di 102 mm per ogni lato, come da normativa

        % Lunghezza del Nose: 
        % Calcolata secondo quanto stabilito dal Sadraey (pag. 429)

        % Lunghezza dei portelloni:
        % stabilita basandosi sulla CS-25.807 Amendment 3 (Book 1, pag. 80)

        % Lunghezza dei bagni:
        % Dal Raymer si ricava il numero di bagni, ogni modulo é lungo un metro

        % Lunghezza della coda:
        % Assunta come la lunghezza del nose, maggiorata di un metro

        % X_dim = [width, length]

        [Business_dim, FC_dim, Economy_dim, Fuselage_Total_Length] = fuselage_size_calculator(Business_geometry, First_class_geometry, Economy_geometry, N_seats);

%% Impennaggi

        % Impennaggio orizzontale

        % Come suggerito dal Sadraey, é stato scelto un profilo NACA 0012
        % per l'impennaggio orizzontale, in quanto simmetrico e non eccessivamente spesso

        % Il coefficiente volumetrico é preso dal Sadraey (pag. 324, tabella 6.4)
        % Il valore ottimale di l é stimato dal Sadraey, prendendo uyn coefficiente k_c = 1 (pag. 324)
        % Il Taper Ratio é compreso tra 0.4 e 0.7, prendiamo 0.6 (é restituito dalla funzione) per ottenere
        % un b_tail di 20 m, paragonabile con i 18 m dell'A350
        % Da considerazioni trigonometriche e, come suggerito dal Sadraey (pag. 339),
        % mantenendo l'angolo di sweep uguale all'ala, si ottiene un determinato impennaggio orizzontale, di seguito plottato

        % Come spiegato dal Sadraey (pag. 340), si parte da un angolo di diedro pari a quello dell'ala, 
        % iterando successivamente (da considerazioni di stabilitá e controllo) fino al valore desiderato

        % L'angolo di calettamento dell'impennaggio orizzontale viene calcolato tramite le formule del Sadraey (pag. 310-311),
        % considerando un volo rettilineo uniforme orizzontale in crociera

        k_c = 1.4;
        Volume_coeff_tail = 1.1;
        C_M0_af = -0.0230;
        Tail_Surface = tail_surface_calculator (k_c, Volume_coeff_tail, MAC, Surface_w, Economy_dim, 'Horizontal');              % [m²]
        alpha_f = 0;
        alpha_h = 1;
        alpha_t = 0;
        i_w = 3.5;
        i_set_tail = set_angle_calculator (alpha_h, alpha_f, i_w);
        
        [TR_tail, Corde_estreme_tail] = wing_profile_plotter(Tail_Surface, Aspect_Ratio, Sweep_Angles, 'Tail');           % Stessi sweep angles dell'ala, stesso gamma
        MAC_tail = mean_aerodynamic_center(Corde_estreme_tail, TR_tail);

        Avg_Mass = average_cruise_mass_calculator (mtow_design, coeff);
        CL_cruise = lift_coefficient_calculator (Avg_Mass, rho_air, Velocity, Surface_w);
        [C_M0_Wing_Body, CL_tail] = horizontal_tail_aero_coefficients_calculator (C_M0_af, alpha_t, Aspect_Ratio, Sweep_Angles, CL_cruise, Volume_coeff_tail);

        % Impennaggio verticale

        % Il coefficiente volumetrico dell'impennaggio verticale é preso 
        % dal Sadraey (pag. 303), nel range suggerito (0.02 % 0.12)
        % Come suggerito dal Sadraey (pag. 322), inizialmente, la distanza tra il centro aerodinamico
        % dell'ala e dell'impennaggio verticale é uguale a quella dell'impennaggio orizzontale

        % L'angolo di calettamento é posto nullo, in quanto il velivolo possiede una configurazione
        % di spinta simmetrica bimotore, a differenza di velivoli con un fan frontale che possono
        % necessitare di un calettamento non nullo, in modo da compensare la coppia di reazione
        % generata dalla singola ventola
        % L'aspect ratio ed il taper ratio sono stati supposti uguali a quelli dell'A350
        % Dal Sadraey (tabella 6.6) si prende uno sweep angle di 35 deg per 
        % l'impennaggio verticale, valore in linea con i velivoli 
        % nella stessa popolazione statistica di quello in studio.
        % Da considerazioni geometriche si ricava poi quello del bordo di fuga.

        i_set_vert = 0;
        Volume_coeff_vert = 0.08;
        Vert_surface = tail_surface_calculator (k_c, Volume_coeff_vert, MAC, Surface_w, Economy_dim, 'Vertical');
        Sweep_Angles_vert = [35, 10.767];

        [TR_vert, Corde_estreme_vert] = wing_profile_plotter(Vert_surface, Aspect_Ratio, Sweep_Angles_vert, 'Vertical');


%% Aerodynamic Feasibility Analysis 

        % L'analisi della portanza alare é molto sensibile alla variazione
        % del taper ratio, in letteratura tipicamente tra 0.3 e 0.4, che risulta
        % essere il range migliore per minimizzare la resistenza indotta, 
        % garantendo quindi la migliore distribuzione di portanza

        % N.B. La funzione della lifting line é implementata a partire da uno script 
        % del Sadraey ed é un modello fisico con forti approssimazioni,
        % per cui é opportuno mettere in evidenza che é adeguato
        % per valutare le prestazioni alari solamente in crociera,
        % non durante manovre critiche o in fase di decollo e atterraggio

        % N.B. Il valore di svergolamento scelto inizialmente era nullo,
        % ma non puó esserlo nell'implementazione numerica del modello, per cui era 
        % inizialmente stato scelto un valore non nullo molto piccolo.
        % Successivamente, si é optato per un valore negativo tra 1 e 3 deg, 
        % in modo da avvicinarsi ad una distribuzione di portanza ellittica, che risulterebbe
        % aerodinamicamente ottimale
        
        Wing_span = wing_span_calculator(Aspect_Ratio, Surface_w);
        CL_wing = lifting_line_plotter (CL_cruise, Corde_estreme, Wing_span, Aspect_Ratio, i_w, Taper_Ratio);
        Lift_wing = lift_calculator (rho_air, Velocity, Surface_w, CL_wing);
        Lift = Lift_wing / 0.8;
        Avg_Weight = Avg_Mass * g;

%% Dimensionamento motore

        % A fronte di una richiesta di spinta (per motore) pari a T = 345.924 kN  viene scelto il seguente motore:

        % General Electric GE90-94B

        % Overall Length: 7.283 [m]
        % Overall WIdth: 3.871 [m]
        % Overall Height: 3.952 [m]
        % Dry Weight: 7892 [kg]
        % Sea Level Static Thrust (Maximum Continuous): 402.920 [kN] 

        % Calcolo delle distanze caratteristiche

        % La Ground Roll Distance é la distanza di rollio in pista, dalla partenza da fermo fino all'istante iniziale del take-off.
        % La Take-Off Distance é la distanza data dalla proiezione della traiettoria di decollo sulla pista, nel punto in cui
        % il velivolo ha superato un ostacolo immaginario di 35 piedi.
        % La BFL (Balanced Field Length) é la lunghezza in cui le distanze percorse in Acc&Go mode e Acc&Stop mode coincidono.
        % La Take Off Run è definita come la distanza percorsa dal punto di inizio decollo al punto medio tra il lift-off e l'ostacolo di 35 piedi.

        Jet_Thrust = 402920;                                                                            % [N]
        N_engines = 2;
        V_min = 57.3803;                                                                                % [m/s]
        friction = 0.05;
        Clearance_Height = 10.67;                                                                       % [m]

        [Ground_Roll_D, Takeoff_D, Takeoff_D_run, Acc_Go_D, Acc_Stop_D] = meaningful_distances_plotter (N_engines, Jet_Thrust, rho_0, V_min, mtow_design, friction, Clearance_Height, Aspect_Ratio, CD_0, Surface_w);
