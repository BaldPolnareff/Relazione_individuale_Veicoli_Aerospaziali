clear all
clc
close all

%% ESERCITAZIONE 1

g = 9.81;                                           % [m / s²]
x = 1e4:1e6;
A = 0.97;
b = -0.06;
y_lett = A .* (x .^ b);
M_TO = [250000 ; 275000 ; 251000 ; 351535 ; 352000 ; 227940 ; 365000 ; 242000 ; 380000 ; 230000 ; 158758 ; 362870 ; 447696 ; 299370 ; 78245 ; 273294 ; 103000]; 
Me = [135500 ; 142400 ; 137000 ; 167829 ; 181400 ; 119950 ; 177755 ; 109400 ; 174000 ; 120600 ; 86069 ; 182480 ; 220128 ; 160530 ; 41413 ; 128808 ; 58300]; 
x1 = M_TO;
y1 = Me ./ M_TO;
figure()
semilogx(x, y_lett, '-', x1(1), y1(1), 'o', x1(2), y1(2), 'o', x1(3), y1(3), 'o', x1(4), y1(4), 'o', x1(5), y1(5), 'o', x1(6), y1(6), 'o', x1(7), y1(7), 'o', x1(8), y1(8), 'o', x1(9), y1(9), 'o', x1(10), y1(10), 'o', x1(11), y1(11), 'o', x1(12), y1(12), 'o', x1(13), y1(13), 'o', x1(14), y1(14), 'o', x1(15), y1(15), 'o', x1(16), y1(16), 'o')
title("Empty mass fraction trend")
xlabel('Takeoff Weight [kg]')
ylabel('Empty weight fraction')
grid on

    %% Miglioramento con dati piú aggiornati, formula del Raymer

n = length(x1);
A1 = [(ones(n, 1)) log(x1)];
c = A1 \ log(y1);
a_new = exp(c(1));
b_new = c(2);
y_new = a_new .* (x .^ b_new);
hold on
semilogx(x, y_new, 'r--')
legend('Statistical trend', 'B-787', 'A350', 'A330 Neo-900', 'B777-300', 'B777 X-900', 'B787-8', 'A340-500', 'A330-3000', 'A340-600', 'A330-200', 'B767-300', 'B747-400', 'B747-8', 'B777-300', 'B737-800', 'MD11', 'Improved trend', 'location', 'SouthWest')
    


    %% Studio preliminare

    % Viene scelta la True Air Speed a 10000 m di quota, in quanto un worse case scenario 
    % rispetto al requisito di Mach 0.85

v = 0.85 * sqrt(1.4 * 287 * (-50 + 273.15));                                        % [m/s]
v = v * 3600;                                                                       % [m/h]
R = 11000 * 10^3;                                                                   % [m]
m_TO = 275000;                                                                      % [kg]           (Valore dell'A350, guess iniziale)
SFC1 = 0.478;                                                                       % [lb/lbh]
Efficienza_max = 20;                                                                % [\]            (Valore ipotizzato)
coeff = [0.97; 0.985; (exp(-R * SFC1 / (v * Efficienza_max))); 0.985; 0.995];
COEFF = prod(coeff);
m_crew = 8 * 85;                                                                    % [kg]           (massa media di 85 kg per persona)

    % Si considera come punto di design un valore di payload max pari a 50 [t]
    % ed un range medio (in crociera) pari a 11000 [km], considerando che, diminuendo in parte
    % il payload, si possono raggiungere range piú lunghi 

m_payload = 50 * 10^3;                                                              % [kg] 
m_to = @(x) x - (m_crew + m_payload) / (1 - 1.06 * (1 - COEFF) - a_new * (x ^ b_new));
m_to_design = fzero(m_to, m_TO);


    %% Valutazione Range: 
    % si valuta come influisce il payload max sulla variazione percentuale del range

range = linspace(10000e+3, 18000e+3, 100);                                          % [km]

for i = 1 : (length(range))
    coeff = [0.97; 0.985; (exp (-range(i) * SFC1 / (v * Efficienza_max))); 0.985; 0.995];
    COEFF = prod(coeff);
    m_to_r = @(x) x - (m_crew + m_payload) / (1 - 1.06 * (1 - COEFF) - a_new * (x^b_new));
    m_to_range(i) = fzero(m_to_r, 275000);
end

figure()
plot(range, m_to_range)
xlabel('range [m]')
ylabel('m_TO [kg]')
grid on

    % Dato che il grafico ha problemi di stabilitá nell'intorno di
    % 15000 [km], si prova ad analizzarne il comportamento 
    % ponendo il payload a 40000 kg

m_payload_2 = 40e+3;                                 % [kg]

for i = 1 : (length(range))
    coeff = [0.97; 0.985; (exp (-range (i) * SFC1 / (v * Efficienza_max))); 0.985; 0.995];
    COEFF = prod(coeff);
    m_to_r = @(x) x - (m_crew + m_payload_2) / (1 - 1.06 * (1 - COEFF) - a_new * (x^b_new));
    m_to_range_2(i) = fzero(m_to_r, 275000);
end

figure()
plot(range, m_to_range_2)
xlabel('range [m]')
ylabel('m_TO [kg]')
grid on

    % Si osserva come un payload minore consenta di raggiungere range maggiori


 
%% VALUTAZIONE PAYLOAD

    % Si considera come range quello di design, pari a 11000 km

R_payload = 11000e+3;                                 % [m]
massa_payload = linspace(20e+3, 80e+3, 100);          % [kg]
coeff = [0.97; 0.985; (exp(-R * SFC1 / (v * Efficienza_max))); 0.985; 0.995];
COEFF = prod(coeff);

for i = 1 : (length(massa_payload))
   m_to_p = @(x) x - (m_crew + massa_payload(i)) / (1 - 1.06 * (1 - COEFF) - a_new * (x^b_new));
   m_to_payload(i) = fzero(m_to_p, 275000);
end

figure()
plot(massa_payload, m_to_payload, 'b')
xlabel('payload [kg]')
ylabel('m_TO [kg]')
grid on


    % Si considera ora come varia il payload in funzione del range, uguale o maggiore di 15000 km

R_payload_2 = 15000e+3;                                % [m]
massa_payload = linspace(20e+3, 80e+3, 100);           % [kg]
coeff = [0.97; 0.985; (exp(-R_payload_2 * SFC1 / (v * Efficienza_max))); 0.985; 0.995];
COEFF = prod(coeff);

for i = 1 : length(massa_payload)
   m_to_p = @(x) x - (m_crew + massa_payload(i)) / (1 - 1.06 * (1 - COEFF) - a_new * (x ^ b_new));
   m_to_payload(i) = fzero(m_to_p, 275000);
end

figure()
plot(massa_payload, m_to_payload, 'b')
xlabel('payload [kg]')
ylabel('m_TO [kg]')
grid on



%% Diagramma Payload-Range

    % Si diagramma il payload in funzione del range del velivolo, considerando come
    % valore di riferimento del fuel quello dell'A350, pari a 141000 L di Avgas.
    % Da questo si calcolano i 4 punti notevoli del grafico PR 

max_fuel_capacity = 141000e-3;                         % [m^3]
rho_avgas = 810;                                       % [kg / m^3]
massa_max_avgas = rho_avgas * max_fuel_capacity;       % [kg]
m_payload_3 = m_to_design * (1 - massa_max_avgas / m_to_design - a_new * (m_to_design ^ b_new)) - m_crew;
Range_3 = -log((1 - massa_max_avgas / m_to_design) / (0.97 * 0.985 * 0.985 * 0.995)) * v * Efficienza_max / SFC1;
m_empty = (a_new * m_to_design ^ b_new) * m_to_design;
m_payload_4 = 0;
m_to_4 = massa_max_avgas + m_empty; 
Range_4 = -log((1 - massa_max_avgas / m_to_4) / (0.97 * 0.985 * 0.985 * 0.995)) * v * Efficienza_max / SFC1;
ypr = [m_payload; m_payload; m_payload_3; m_payload_4];
xpr = [0; R; Range_3; Range_4] * 10 ^(-3);
figure()
hold on
plot(xpr, ypr, 'r')
xlabel('range [km]')
ylabel('payload [kg]')
title('Payload - Range')
grid on
hold on

    % Viene considerata, per diagrammare la famiglia di velivoli, una
    % variazione del MTOW pari a 0.88 * m_to_design, per il velivolo con piú
    % alto range, con un payload piú piccolo di 5e+3 kg 
    % ed una variazione per il velivolo con piú alto 
    % payload di MTOW pari a 1.077 * m_to_design. 
    % Con un payload maggiore di 5e+3 kg
    % il procedimento seguito é sempre lo stesso. 
    % Viene fatto variare anche il range del punto 2 ponendolo ± 500000 m

m_payload_1 = 45e+3;
m_to_design1 = m_to_design * 0.88;
R_1 = 11500000;
m_payload_3 = m_to_design1 * (1 - massa_max_avgas / m_to_design1 - a_new * (m_to_design1 ^ b_new)) - m_crew;
Range_3 = -log((1 - massa_max_avgas / m_to_design1) / (0.97 *0.985 *0.985 *0.995)) * v * Efficienza_max / SFC1;
m_empty = (a_new * m_to_design1 ^ b_new) * m_to_design1;
m_payload_4 = 0;
m_to_4 = massa_max_avgas + m_empty; 
Range_4 = -log((1 - massa_max_avgas / m_to_4) / (0.97 * 0.985 * 0.985 * 0.995)) * v * Efficienza_max / SFC1;
ypr = [m_payload_1; m_payload_1; m_payload_3; m_payload_4];
xpr = [0; R_1; Range_3; Range_4] * 10 ^(-3);
plot(xpr, ypr, 'b')
xlabel('range [km]')
ylabel('payload [kg]')
title('Payload - Range')
grid on

m_payload_2 = 55e+3;
R_2 = 10500000;
m_to_design_2 = m_to_design * 1.077;
m_payload_3 = m_to_design_2 * (1 - massa_max_avgas / m_to_design_2 - a_new * (m_to_design_2 ^ b_new)) - m_crew;
Range_3 = -log((1 - massa_max_avgas / m_to_design_2) / (0.97 * 0.985 * 0.985 * 0.995)) * v * Efficienza_max / SFC1;
m_empty = (a_new * m_to_design_2 ^ b_new) * m_to_design_2;
m_payload_4 = 0;
m_to_4 = massa_max_avgas + m_empty; 
Range_4 = -log((1 - massa_max_avgas / m_to_4) / (0.97 * 0.985 * 0.985 * 0.995)) * v * Efficienza_max / SFC1;
ypr = [m_payload_2; m_payload_2; m_payload_3; m_payload_4];
xpr = [0; R_2; Range_3; Range_4] * 10^(-3);
plot(xpr, ypr, 'g')
xlabel('range [km]')
ylabel('payload [kg]')
title('Payload - Range')
grid on
axis([0 2.5e+4 0 6.5e+4])
legend('M_TOW = 304685 [kg]','M_TOW = 268123 [kg]','M_TOW = 328146 [kg]','Location','NorthEast')






%% MATCHING CHART


rho_0 = 1.225;                                      % [kg / m^3]
T0 = 288;                                           % [K]
z = 12000;                                          % [m]
h = -0.0065;                                        % [K / m]
rho = rho_0 * ((T0 + h * z) / T0) ^4.2561;          % [kg / m^3]

    % Passo 1 -> Velocitá di stallo %

    % Il Raymer stabilisce di prendere la velocitá di approccio da velivoli simili,
    % in quanto non presente nelle normative vigenti: 
    % vengono quindi presi due valori, simili a quello dell'A350-900, pari a 140 kts e
    % a quello dell'A350-1000, pari a 147 kts

v_approach = 145 * 0.514444;            % [m / s]
v_stallo = v_approach / 1.3;            % [m / s]


    % Una volta consumato il fuel necessario (per il range di 11000 km), 
    % la massa al landing risulta essere circa 214000 kg, valore in linea
    % con i velivoli di riferimento (circa tra 205000 kg e 236000 kg)

    % Secondo quando riportato dal Raymer, valori tipici per il 
    % coefficiente di portanza massimo di un jet transport sono
    % all'incirca tra 2.2 e 3.2.
    % É stato quindi scelto un valore pari a 2.8

    % Prendendo ancora una volta spunto dal Raymer, i calcoli sono svolti
    % con riferimento alla densitá al livello del mare. 

cl_max = 2.8;

W_S = 0.5 * rho_0 * (v_stallo ^2) * cl_max / g;

figure()
hold on
line("xdata", [W_S, W_S], "ydata", [-1000, 1000], "linewidth", 2)
xlabel('W/S [kg / m²]')
ylabel('T/W')
title('Matching Chart') 
grid on






    % Passo 2 -> Requisiti di cruise 

    % N.B. La W/S riportata nei grafici di riferimento esprime la forza
    % in kilogrammi peso, per cui sono state fatte le dovute correzioni dimensionali, 
    % riportando le forze in newton, in modo da ottenere (correttamente)
    % un rapporto spinta-peso adimensionale.


v = 0.85 * sqrt(1.4 * 287 * (-50 + 273.15));        % [m / s]
v_max = 1.25 * v;                                   % [m / s]
C_D_0 = 0.015;
AR = 9;
e = 0.921;
K = 1 / (pi * e * AR);
sigma = rho / rho_0;

x2 = linspace(0, 700);                              % [kg / m²]
y2 = @(x) rho_0 * C_D_0 * (v_max ^2) ./ (2 * g * x) + 2 * g * K * x ./ (rho * sigma * (v_max ^2));

plot(x2, y2(x2), 'LineWidth', 1.5)
y2(W_S);

    % Passo 3 -> Take Off distance %

    % É stata considerata una ground roll distance di 888 m (circa 3000 ft).
    % La "over 50 ft" distance considerata é maggiore di 1.7 volte.
    % Tale valore é quello a cui corrisponde un TOP (Take Off Parameter)
    % di circa 230 lbs/ft².
    % Inoltre, come suggerito dal Roskam (pag. 107), il coefficiente di 
    % portanza al take off é pari a quello massimo, diviso per un fattore di 1.21


TOP = 230 * 0.45 / (0.3048 ^ 2);                     % [kg / m²]
cl_max_to = 2.5;
cl_to = cl_max_to / 1.21;                

y3 = @(x) (x) ./ (TOP * cl_to);

plot(x2, y3(x2), 'LineWidth', 1.5)
axis([0 700 0 1])

Efficienza_max = 20;



    % Passo 4 %

C_D_0_to = 0.075 + 0.02;
ROC = 3000 * 0.00508;                                % [m / s]
y4 = @(x) ROC ./ sqrt(2 * g * x ./ (rho_0 * sqrt(C_D_0_to / K))) + 1 / Efficienza_max;
plot(x2, y4(x2), 'LineWidth', 1.5)


    % Passo 5 %

line("xdata", [-1000, 1000], "ydata", [1 / (sigma * Efficienza_max), 1 / (sigma * Efficienza_max)], "linewidth", 3)




%% Punto di design

plot(W_S, y3(W_S), '*r', 'LineWidth', 3)
legend('STALLO', 'CRUISE', 'TOP', 'ROC', 'CEILING', 'Design Point', 'Location', 'NorthEast')

    % Il punto di design scelto é y3 che, come visibile nel matching chart,
    % é il punto in grado di minimizzare la superficie alare, a discapito di un
    % leggero aumento della spinta necessaria.
    % Tale scelta é motivata dalla necessitá di ridurre la superficie alare,
    % inizialmente troppo alta, dato che il motore scelto ha una spinta
    % disponibile piú che sufficiente a coprire la spinta necessaria,
    % anche con un ampio margine.









%% PROGETTO ALA


    % Spinta di progetto

T = y3(W_S) * m_to_design * g;                                                  % [N] 
S = m_to_design / W_S;                                                          % [m²]

    % N.B. Si fa riferimento ad una massa media tra inizio e fine crociera
    
m_in_c = coeff(1) * coeff(2) * m_to_design;
m_fin_c = coeff(1) * coeff(2) * coeff(3) * m_to_design;

m_avg_tmp = (m_in_c + m_fin_c) / 2;

m_avg = m_avg_tmp;
C_L_cruise = 2 * m_avg * g / (rho * v ^2 * S);
k_w = 0.95;                                                                     % Tipicamente 95% della spinta verticale
C_L_cruise_w = C_L_cruise / k_w; 
k_a = 0.9;
C_L_i = C_L_cruise_w / k_a;
CL_MAX = 2 * m_to_design * g / (rho_0 * (v_stallo * 1.1) ^2 * S);
CL_MAX_W = CL_MAX / k_w;
CL_MAX_gross = CL_MAX_W / k_a;

    % Secondo indicazione del Sadraey, é stata ipotizzata la presenza di 
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
    % Una volta calcolata la corda media aerodinamica (vedi sezioni successive), si potrá
    % verificare il Reynolds ottenuto e validare (o scartare) i dati presi dal sito.
    

 %% Dimensionamento geometrico dell'ala 

b = sqrt(AR * S);                                                                % [m]
S_w_ang = 29;                                                                    % [deg]
S_w_ang_te = 18;                                                                 % [deg]
Gamma = 5;                                                                       % [deg]
b_med = b / 2;                                                                   % [m]
c_root_wing = (S / (b_med) + b_med * (tand(S_w_ang) - tand(S_w_ang_te))) / 2;      % [m]
c_tip_wing = c_root_wing - b_med * (tand(S_w_ang) - tand(S_w_ang_te));              % [m]
TR = c_tip_wing / c_root_wing;                                                      % (Tapering Ratio)
y1 = c_root_wing - b_med * tand(S_w_ang);
y2 = y1 - c_tip_wing;
figure()
plot([0 b_med], [0 y2])
hold on
plot([0 b_med], [c_root_wing y1])
plot([0 0], [0 c_root_wing])
plot([b_med b_med], [y2 y1])
grid on
xlabel('asse y')
ylabel('asse x')






















%% ES4

    % I parametri sono stati scelti seguendo i suggerimenti del Sadraey (pag. 376)

n_business = 26;                                                                % Posti in business
n_fc = 48;                                                                      % Posti in first class
n_ec = 315;                                                                     % Posti in economy (9 per fila)
W_pax_ec = n_ec * (82 + 23 + 14);                                               % Massa raccomandata per ogni passeggero in economy
W_pax_oth = (n_fc + n_business) * (82 + 32 + 14);                               % Massa raccomandata per ogni passeggero in first class e business


    % Configurazione scelta

    % Business:     1-2-1
    % First Class:  2-2-2 
    % Economy:      3-3-3

wa_bn = 70e-2;                                                                  % [m] (larghezza corridoio)
wa_fc = 60e-2;                                                                  % [m]
wa_ec = 50e-2;                                                                  % [m]
ps_ec = 72e-2;                                                                  % [m] (lunghezza tra sedili)
ps_fc = 92e-2;                                                                  % [m]
ps_bn = 104e-2;                                                                 % [m]
ws_ec = 46e-2;                                                                  % [m] (larghezza sedili)
ws_fc = 60e-2;                                                                  % [m]
ws_bn = 75e-2;                                                                  % [m]





    % BUSINESS

w_bn = 4 * ws_bn + 2 * wa_bn;                                                           % [m]
l_bn = n_business / 4 * ps_bn;                                                          % [m]

    % FIRST CLASS

w_fc = 6 * ws_fc + 2 * wa_fc;                                                           % [m]
l_fc = n_fc / 6 * ps_fc;                                                                % [m]

    % ECONOMY

w_ec = 9 * ws_ec + 2 * wa_ec;                                                           % [m]
l_ec = n_ec / 9 * ps_ec;                                                                % [m]


    % Vengono aggiunti 102 mm di spessore per ogni lato, come stabilito da normativa


   
    % Lunghezza del Nose 
    % Calcolata secondo quanto stabilito dal Sadraey (pag. 429)

l_nose = 1.75 * (w_ec + 2 * 0.102);                                                     % [m]

    % Lunghezza dei portelloni:
    % stabilita basandosi sulla CS-25.807 Amendment 3 (Book 1, pag. 80)

l_portelloni = 4 * 1.07;                                                                % [m]

    % Lunghezza dei bagni:
    % Dal Raymer si ricava il numero di bagni, ogni modulo é lungo un metro

l_bath = 4;                                                                             % [m]

    % Lunghezza della coda

l_coda = l_nose + 1;                                                                    % [m]


    % lunghezza totale

l_tot = l_bn + l_fc + l_ec +l_nose + l_portelloni + l_bath + l_coda;                    % [m]













%% ESERCITAZIONE 5  

    % Come suggerito dal Sadraey, é stato scelto un profilo NACA 0012
    % per l'impennaggio orizzontale, in quanto simmetrico e non eccessivamente spesso


    % La cordia media aerodinamica viene calcolata con una formula approssimata,
    % fornita dal Raymer (pag. 104) ed in seguito usata anche per gli impennaggi

MAC = 2 / 3 * c_root_wing * (1 + TR + TR ^ 2) / (1 + TR);

    % Il coefficiente volumetrico é preso dal Sadraey (pag. 324, tabella 6.4)

V_h = 1.1;

    % Il valore ottimale di l é stimato dal Sadraey, prendendo uyn coefficiente k_c = 1 (pag. 324)

k_c = 1.4;
l_opt = k_c * sqrt(4 * MAC * S * V_h / (pi * w_ec));
S_tail = V_h * MAC * S / l_opt;


    % Il Taper Ratio é compreso tra 0.4 e 0.7, prendiamo 0.6 per ottenere
    % un b_tail di 20 m, paragonabile con i 18 m dell'A350

TR_tail = 0.6;

    % Da considerazioni trigonometriche e, come suggerito dal Sadraey (pag. 339),
    % mantenendo l'angolo di sweep uguale all'ala, si ottengono i seguenti risultati:

Sweep_tail = 29;
Sweep_te_tail = 18;
b_tail = 2 * sqrt(S_tail * (1 - TR_tail) / ((1 + TR_tail) * (tand(Sweep_tail) - tand(Sweep_te_tail))));
c_root_tail = (b_tail / 2) * (tand(Sweep_tail) - tand(Sweep_te_tail)) / (1 - TR_tail);
c_tip_tail = TR_tail * c_root_tail;
AR_tail = b_tail ^ 2 / S_tail;

    % Controllo sul Reynolds

delta_CL_HLD = 1.4;                                                                % prendiamo 1 di Fowler Flap e 0.3 di Leading Edge slat dal Sadrey ( prendiamo i valori più bassi per essere conservativi)
Cl_max = CL_MAX_gross - delta_CL_HLD;
dyn_viscosity = 1.46e-5;
Reynolds = rho * v * MAC / dyn_viscosity; 

    % Come si puó vedere, i valori di AirfoilTools sono stati validati.











    % Aerodinamica dell'impennaggio orizzontale
C_M0_af = -0.0230;                                                                  % (Da airofoiltools NACA1412)
alpha_t = 0;                                                                        % (Angolo di twist)
C_M0_wf = C_M0_af * (AR * cosd(S_w_ang) ^ 2) / (AR + 2 * cosd(S_w_ang)) + 0.01 * alpha_t;
CL_H = (C_M0_wf + C_L_cruise * (0.2)) / V_h;                                        % Valore medio del delta h preso dal Sadraey
l_h = l_opt - 0.2 * MAC;
MAC_t = 2 / 3 * c_root_tail * (1 + TR_tail + TR_tail ^ 2) / (1 + TR_tail);
Gamma_t = 5;                                                                        % [deg] 

    % Come spiegato dal Sadraey (pag. 340), si parte da un angolo di diedro pari a quello dell'ala, 
    % iterando successivamente (da considerazioni di stabilitá e controllo) fino al valore desiderato

    % L'angolo di calettamento dell'impennaggio orizzontale viene calcolato tramite le formule del Sadraey (pag. 310-311),
    % considerando un volo rettilineo uniforme orizzontale in crociera

alfa_h = 1;                                                                         % [deg]
alfa_f = 0;                                                                         % [deg]
i_w = 3.5;                                                                          % [deg]
alfa_w = i_w + alfa_f;                                                              % [deg]
eps_0 = 1;                                        
deps_dalfa = 0.3;
eps = eps_0 + deps_dalfa * alfa_w;
i_h = alfa_h - alfa_f + eps;                                                        % [deg]



%% Dimensionamento dell'impennaggio verticale

    % Il coefficiente volumetrico dell'impennaggio verticale é preso 
    % dal Sadraey (pag. 303), nel range suggerito (0.02 % 0.12)

V_vtail = 0.08;

    % Come suggerito dal Sadraey (pag. 322), inizialmente, la distanza tra il centro aerodinamico
    % dell'ala e dell'impennaggio verticale é uguale a quella dell'impennaggio orizzontale

l_vtail = l_opt;
S_vtail = V_vtail * (S * b) / l_vtail;

    % L'angolo di calettamento é posto nullo, in quanto il velivolo possiede una configurazione
    % di spinta simmetrica bimotore, a differenza di velivoli con un fan frontale che possono
    % necessitare di un calettamento non nullo, in modo da compensare la coppia di reazione
    % generata dalla singola ventola

i_vtail = 0;

    % L'aspect ratio ed il taper ratio sono stati supposti uguali a quelli dell'A350
 
TR_vtail = 0.3902;
AR_vtail = 1.72;
b_vtail = sqrt(AR_vtail * S_vtail);
c_root_vtail = 2 * S_vtail / ((1 + TR_vtail) * b_vtail);
c_tip_vtail = TR_vtail * c_root_vtail;

    % Dal Sadraey (tabella 6.6) si prende uno sweep angle di 35 deg per 
    % l'impennaggio verticale, valore in linea con i velivoli 
    % nella stessa popolazione statistica di quello in studio.
    % Da considerazioni geometriche si ricava poi quello del bordo di fuga.

Sweep_vtail = 35;                                                                                       % [deg]
Sweep_vtail_te = (180 / pi) * atan((c_tip_vtail - c_root_vtail + b_vtail * tand(35)) / b_vtail);        % [deg]

    % Si riportano i valori relativi al vertical tail dell'A350:
    % c_root_rudder = 7.79;                                                                             % [m]
    % c_tip_rudder = 3.04;                                                                              % [m]
    % TR_rudder = c_tip_rudder / c_root_rudder;
    % b_rudder = 9.42;                                                                                  % [m]
    % S_rudder = (c_root_rudder + c_tip_rudder) * b_rudder / 2;                                         % [m²]




%% Feasibility analysis -> Aerodynamic

    % N.B. L'analisi della portanza d'ala é molto sensibile alla variazione
    % di taper ratio, in letteratura tipicamente tra 0.3 e 0.4, che risulta
    % essere il range migliore per minimizzare la resistenza indotta, 
    % garantendo quindi la migliore distribuzione di portanza

    % N.B. Il codice seguente é preso dal Sadraey ed é un modello fisico
    % con forti approssimazioni, per cui é opportuno mettere in evidenza
    % che é adeguato per valutare le prestazioni alari solamente
    % in crociera, non durante manovre critiche o in fase di 
    % decollo e atterraggio

N_segmenti = 200;                                                                                       % [number of segments - 1]
alpha_twist = -2;                                                                                       % [deg] (Angolo di twist)
a_2d = 6.3;                                                                                             % [1 / rad] (lift curve slope)
alpha_0 = -1.5;                                                                                         % [deg] (flap down zero-lift angle of attack)
theta = pi / (2 * N_segmenti) : pi / (2 * N_segmenti) : pi/2;
alpha = i_w + alpha_twist : -alpha_twist / (N_segmenti - 1) : i_w;
    % Angolo d'attacco del segmento
z = (b / 2) * cos(theta);
MAC_segmento = c_root_wing * (1 - (1 - TR) * cos(theta));                                                 % MAC ad ogni segmento
mu = MAC_segmento * a_2d / (4 * b);
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
CL = 4 * b * sum2 ./ MAC_segmento;
CL_1 = [0, CL];
y_s = [b / 2, z];
figure()
plot(y_s,CL_1,'-o')
grid on
xlabel('semiala [m]')
ylabel('Cl - Lift Coefficient')
title('Andamento Cl - Teoria della Linea Portante')
CL_wing = pi * AR * AA(1);                                                                              % Coefficiente di portanza di tutta l'ala

    % Il contributo di portanza dato dall'ala é circa l'85% del totale

if CL_wing >= 0.85 * m_avg * g * 2 / (rho * v ^ 2 * S)
    display('Corretto')
else 
    display('Errore: Rivedere i valori inseriti')
end

    % N.B. Il valore di svergolamento scelto inizialmente era nullo,
    % ma non puó esserlo nella concezione numerica del modello, per cui era inizialmente 
    % stato scelto un valore non nullo molto piccolo.
    % Successivamente, si é optato per un valore negativo 1 e 3 deg, 
    % in modo da avvicinarsi ad una distribuzione di portanza ellittica, che risulterebbe
    % aerodinamicamente ottimale

    % Portanza e forza peso

L = 0.5 * rho * v ^ 2 * S * CL_wing / 0.8;
W = m_avg * g;





%% Dimensionamento motore

    % A fronte di una richiesta di spinta pari a T = 345.924 kN
    % viene scelto il seguente motore:


    %% Rolls Roy_ce Trent XWB-97

    % Lunghezza = 5.812 m
    % Massa = 7550 kg
    % Take-off Thrust = 431 kN
    % Diametro Fan = 3.00 m

W = m_to_design * g;                                                                    % [N]
CL_TO = 2;                                                                              % Coefficiente di portanza in take off
    % T_12 = 379000;                                                                    % [N] (Trent XWB 84)
T_12 = 360400;                                                                          % [N] (TRENT 1000-R3 (2017))
N_engine = 2;                                                                           % Numero motori
T = T_12 * N_engine;                                                                    % Total thrust
friction = 0.05;                                                                        % Coefficiente di attrito
h_0 = 10.67;                                                                            % Distanza ostacolo in metri
k = 1 / (pi * e * AR);
C_lift = friction / (2 * k);
n_x = 0.3;                                                                              % Fattore di carico lungo x

    % Decollo

V_min = 57.3803;
V_R = 1.2 * V_min;                                                                      % [m/s] Velocità di rotazione
A_gr = g * (T / W - friction);
B = (rho * g) / (2 * W / S) * (C_D_0 - (friction ^ 2) / (4 * k));
d_roll = 1 / (2 * B) * (log(A_gr) - log(A_gr - B * (V_R ^ 2)));
d_liftoff = 2 * V_R;
n_z = (0.734 * V_R ^ 2) / (V_min ^ 2); 
cosgamma = 1 - (h_0) / R;
gamma = acos(cosgamma); 
d_H = R * sin(gamma);
d_takeoff = d_roll + d_liftoff + d_H;

    % BFL

V_fail = linspace(0, V_R, 1000);
d_aeo = 1 / (2 * B) * (log(A_gr) - log(A_gr - B * (V_fail .^ 2)));

    % Nel tratto successivo vi é una variazione di A

A_oei = g * (0.5 * T / W - friction);
d_oei = 1 / (2 * B) * (log(A_oei - B * (V_fail .^ 2)) - log(A_oei - B * (V_R ^ 2)));
d_AAG = d_aeo + d_oei + d_H + d_liftoff;
d_react = 2 * V_fail;
d_stop = (V_fail .^ 2) / (2 * g * n_x);
d_AAS = d_aeo + d_react + d_stop;
figure()
plot(V_fail, d_AAG, 'r', 'linewidth', 1.5)
hold on
plot(V_fail, d_AAS, 'b', 'linewidth', 1.5)
grid on
title('Calcolo della distanza di BFL')
xlabel('Velocità della piantata motore V_F_A_I_L[m/s]')
ylabel('Distanza [m]')
legend('Acceleration & GO','Acceleration & STOP')

    % L'ultimo passo da fare è ricavare la Takeoff Run, 
    % definita dalla normativa CS-25.113C (Emendamento 2)


d_TO_run = (d_roll + d_liftoff + d_H / 2);

    % Vengono plottate la Takeoff Distance e la Goll

x_c = d_roll + d_liftoff;
y_c = R;
xx = linspace(x_c, d_takeoff);
yy = y_c - sqrt(R ^ 2 - (xx - x_c) .^ 2);
xx_1 = linspace(0, x_c);
yy_1 = zeros(1, length(xx_1));

figure()
plot([0 d_TO_run], [0 0], 'LineWidth', 1.5)
grid on
xlabel('Distanza [m]')
ylabel('Quota [m]')
title('Take Off Run')


figure()
plot([0 d_takeoff], [0 0], 'LineWidth', 1.5)
grid on
xlabel('Distanza [m]')
ylabel('Quota [m]')
title('Take Off Distance')

figure()
plot([xx_1, xx], [yy_1, yy], 'LineWidth', 1.5)
grid on
xlabel('Distanza [m]')
ylabel('Quota [m]')
title('Traiettoria effettiva del velivolo')

figure()
hold on
plot([xx_1, xx], [yy_1, yy], 'LineWidth', 1.5)
plot([0 d_takeoff], [0 0], 'LineWidth', 1.5)
plot([0 d_TO_run], [0 0], 'LineWidth', 1)
grid on
xlabel('Distanza [m]')
ylabel('Quota [m]')
title('Take Off Distance')
legend('Traiettoria effettiva del velivolo','Take off Distance','Take Off Run','Location','northeast')


