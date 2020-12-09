clear all
clc
close all

% N.B. Le variabili esplicitamente temporanee (TMP_i) sono da considerarsi come tali e fini solo alla forma

%% ESERCITAZIONE 1  

    % Variabili fisse e variabili di statistical trend 

        global g = 9.81;                                               % [m/s²]
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
        % Per problemi di stabilitá del plot nell'intorno di 15000 km, é stato poi scelto un payload di 40 tonnellate.

        velocity = 3600 * 0.85 * sqrt(1.4 * 287 * (-50 + 273.15));      % [m/h]
        Range = 11000 * 10 ^ 3;                                         % [m]
        SFC_1 = 0.478;                                                  % [lb/lbh]
        m_TO_guess = 275000;                                            % [kg]           (Valore del velivolo target (A350), guess iniziale)
        m_payload = 40e+3;                                              % [kg]
        N_crew = 8;
        Efficienza_max = 20;

        [MTOW] = range_plotter (velocity, Range, N_crew, SFC_1, Efficienza_max, m_payload, m_TO_guess, a_new, b_new);
        mtow_design = MTOW(1);
        mtow_range = MTOW(2);

    % Valutazione Payload 

        % Si considera il range di design, pari a 11000 km

        