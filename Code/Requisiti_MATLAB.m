clear all
clc
close all
x = 1e4:1e6;
A=0.97;
b=-0.06;
y_lett = A.*(x.^b);
MTO = [250000 ; 275000 ; 251000 ; 351535 ; 352000 ; 227940 ; 365000 ; 242000 ;380000 ; 230000 ;	158758 ; 362870 ; 447696 ; 299370 ; 78245 ; 273294  ; 103000]; 
Me = [135500 ; 142400 ; 137000 ; 167829 ; 181400 ; 119950 ; 177755 ; 109400 ; 174000 ; 120600 ; 86069 ; 182480 ; 220128 ; 160530 ; 41413 ; 128808 ; 58300]; 
x1 = MTO;
y1 = Me./MTO;
figure()
semilogx(x,y_lett,'-',x1(1),y1(1),'o',x1(2),y1(2),'o',x1(3),y1(3),'o',x1(4),y1(4),'o',x1(5),y1(5),'o',x1(6),y1(6),'o',x1(7),y1(7),'o',x1(8),y1(8),'o',x1(9),y1(9),'o',x1(10),y1(10),'o',x1(11),y1(11),'o',x1(12),y1(12),'o',x1(13),y1(13),'o',x1(14),y1(14),'o',x1(15),y1(15),'o',x1(16),y1(16),'o')
title("Empty mass fraction trend")
xlabel('Takeoff Weight [kg]')
ylabel('Empty weight fraction')
grid on

%% Miglioramento con dati pi� aggiornati formula del Rymer
n=length(x1);
A1=[ones(n,1) log(x1)];
c=A1\log(y1);
a_new=exp(c(1));
b_new=c(2);
y_new= a_new.*(x.^b_new);
hold on
semilogx(x,y_new,'r--')  
legend('Trend equation','B-787','A350','A330 Neo-900','B777-300','B777 X-900','B787-8','A340-500','Improved trend','location','SouthWest')

% %% Grafico Efficienza
% A=[10.03,9.49,11,9.82,9.96];
% Sref=[377,443,465,436.8,516.7];
% figure(2)
% LD=15.5*sqrt(A/6.2);
% plot(A(1)/6.2,LD(1),'*',A(2)/6.2,LD(2),'*',A(3)/6.2,LD(3),'*',A(4)/6.2,LD(4),'*',A(5)/6.2,LD(5),'*')
% title("Max Efficiency")
% xlabel('Wetted Aspect Ratio [=A/(Swet/Sref)]')
% ylabel('L/D max')
% legend('B-787','A350','A330 Neo-900','B777-300','B777 X-900','location','SouthEast')
% grid on
% 
% %% Grafico SFC
% SFC=[0.506,0.478,0.506,0.56,0.545];
% M_cruise=[0.85,0.85,0.81,0.84,0.84];
% figure(3)
% plot(M_cruise(1),SFC(1),'*',M_cruise(2),SFC(2),'*',M_cruise(3),SFC(3),'*',M_cruise(4),SFC(4),'*',M_cruise(5),SFC(5),'*')
% legend('B-787','A350','A330 Neo-900','B777-300','B777 X-900','location','best')
% title("Equivalent Jet SFC")
% xlabel('Mach')
% ylabel('SFC [lb/lbh]')
% legend('B-787','A350','A330 Neo-900','B777-300','B777 X-900','location','best')
% grid on
% xlim([0.8,0.9]);
% ylim([0.4,0.7]);

%% Studio preliminare
%velocit� usata TAS, quota 10000 metri perch� condizione peggiore del
%requisito, Mach 0.85

v=0.85*sqrt(1.4*287*(-50+273.15)); %m/s
v=v*3600; %m/h
R=11000*10^3; %m
mTO=275000; %[kg] A350, lo uso come guess
SFC1=0.478;%lb/lbh
Efficienza_max=20; %Valore ipotizzato
coeff=[0.97;0.985;exp(-R*SFC1/(v*Efficienza_max));0.985;0.995];
COEFF=prod(coeff);
mcrew=8*85; %[kg]  massa media  di 85 kg 

% si considera come punto di design un valore di payload max pari a 55 [t]
% e un range medio pari a 11000 [km], considerando che diminuendo in parte
% il payload si possono raggiungere range pi� lunghi (vedi dopo)

mpayload=55*10^3; %[kg] 
mto=@(x) x-(mcrew+mpayload)/(1-1.06*(1-COEFF)-a_new*(x^b_new));
mto_design=fzero(mto,mTO);
display(mto_design)


%% Valutazione Range 
% si valuta come influisce il payload max sulla variazione percentuale del range

range=linspace(10000e+3,18000e+3,100);  %[km]
for i=1:length(range)
    coeff=[0.97;0.985;exp(-range(i)*SFC1/(v*Efficienza_max));0.985;0.995];
    COEFF=prod(coeff);
    mto_r=@(x) x-(mcrew+mpayload)/(1-1.06*(1-COEFF)-a_new*(x^b_new));
    mto_range(i)=fzero(mto_r,275000);
end
figure()
plot(range,mto_range)
xlabel('range [m]')
ylabel('mTO [kg]')
grid on

% dato che il grafico da problemi di stabilit� intorno ad un valore di
% 15000 [km] si prova a vedere come si comporta diminuendo il payload e
% ponendolo a 40 [t]

mpayload2=40e+3;  % [kg]
for i=1:length(range)
    coeff=[0.97;0.985;exp(-range(i)*SFC1/(v*Efficienza_max));0.985;0.995];
    COEFF=prod(coeff);
    mto_r=@(x) x-(mcrew+mpayload2)/(1-1.06*(1-COEFF)-a_new*(x^b_new));
    mto_range2(i)=fzero(mto_r,275000);
end
figure()
plot(range,mto_range2)
xlabel('range [m]')
ylabel('mTO [kg]')
grid on

% si osserva come giustamente con un payload minore si raggiungano range
% maggiori 

%% Valutazione Payload
% si considera come valore di range quello di design e pari a 11000 [km]
R_payload=11000e+3;  %[m]
massa_payload=linspace(20e+3,80e+3,100); %[kg]
coeff=[0.97;0.985;exp(-R*SFC1/(v*Efficienza_max));0.985;0.995];
COEFF=prod(coeff);
for i=1:length(massa_payload)
   mto_p=@(x) x-(mcrew+massa_payload(i))/(1-1.06*(1-COEFF)-a_new*(x^b_new));
   mto_payload(i)=fzero(mto_p,275000);
end
figure()
plot(massa_payload,mto_payload,'b')
xlabel('payload [kg]')
ylabel('mTO [kg]')
grid on


% si considera ora come varia il payload in funzione di un valore di range
% maggiore e pari a 15000 [km]

R_payload2=15000e+3;  %[m]
massa_payload=linspace(20e+3,80e+3,100); %[kg]
coeff=[0.97;0.985;exp(-R_payload2*SFC1/(v*Efficienza_max));0.985;0.995];
COEFF=prod(coeff);
for i=1:length(massa_payload)
   mto_p=@(x) x-(mcrew+massa_payload(i))/(1-1.06*(1-COEFF)-a_new*(x^b_new));
   mto_payload(i)=fzero(mto_p,275000);
end
figure()
plot(massa_payload,mto_payload,'b')
xlabel('payload [kg]')
ylabel('mTO [kg]')
grid on


%% Diagramma Payload-Range
% si calcola il diagramma payload range del velivolo considerando come
% valore di riferimento del fuel quello dell'a350 pari a 141000 l di avgas.
% Da questo si calcolano i 4 punti notevoli del grafico PR 
max_fuel_capacity = 141000e-3; %[m^3]
rho_avgas = 810;               %[kg/m^3]
massa_max_avgas = rho_avgas*max_fuel_capacity;  %[kg]
mpayload3 = mto_design*(1-massa_max_avgas/mto_design-a_new*(mto_design^b_new))-mcrew;
range3= -log((1-massa_max_avgas/mto_design)/(0.97*0.985*0.985*0.995))*v*Efficienza_max/SFC1;
m_empty = (a_new*mto_design^b_new)*mto_design;
mpayload4 = 0;
mto4 = massa_max_avgas + m_empty; 
range4 = -log((1-massa_max_avgas/mto4)/(0.97*0.985*0.985*0.995))*v*Efficienza_max/SFC1;
ypr = [mpayload ; mpayload ; mpayload3 ; mpayload4];
xpr = [0 ; R ; range3 ; range4]*10^(-3);
figure()
plot(xpr,ypr,'r')
xlabel('range [km]')
ylabel('payload [kg]')
title('Grafico Payload - Range')
grid on

%% AGGIUNGERE LA PARTE RELATIVA ALLA FAMIGLIA DI VELIVOLI!!!!! %%

%% Matching Chart
clear all 
clc
close all

rho0 = 1.225;                         % [kg/m^3]
T0 = 288;                             % [K]
z = 12000;                            % [m]
h = -0.0065;                          % [K/m]
rho = rho0*((T0+h*z)/T0)^4.2561;      % [kg/m^3]

% Passo 1 -> Velocit� di stallo %

% il Raymer dice di prendere la velocit� di approccio dai velivoli simili
% in quanto non � presente nelle normative vigenti, allora noi abbiamo
% preso un valore simile a quello dell'a350-900 pari a V_approccio = 140 kts e
% quello dell'a350-1000 = 147 kts

v_approach = 140*0.514444;            % [m/s]
v_stallo = v_approach/1.3;            % [m/s]

% Si � controllato come la massa che avrebbe il velivolo al landing,
% consumando il fuel necessario per il range di riferimento (R=11000 [km])
% veniva 214 [t] che � inlinea con i velivoli di riferimento che erano tra
% le 205 e 236 [ton]

% Dal Raymer si � ricavato che valori tipici per il cl_max di un jet
% transport stanno tra 2.2 e 3.2 quindi noi scegliamo un valore pari a 3

% si utilizza dal Raymer il valore della rho al take off ovvero rho0 =
% 1.225 kg/m^3

cl_max = 3;

W_S = 0.5*rho0*(v_stallo^2)*cl_max/9.81;

figure()
hold on
xline(W_S,'b','Linewidth',1.5)
xlabel('W/S [kg/m^2]')
ylabel('T/W')
title('Matching Chart') 
grid on

% Passo 2 -> Requisiti di cruise 

% Attenzione che in questo caso per rendere adimensionale la T/W si
% moltiplica la W/S che sul grafico � considerata in kg/m^2 per la
% accelerazione di gravit� considerata come 9.81

v = 0.85*sqrt(1.4*287*(-50+273.15));
v_max = 1.25*v;
cd0 = 0.015;
AR = 12;
e = 0.921;
K = 1/(pi*e*AR);
sigma = rho/rho0;

x2 = linspace(0,700);            % [kg/m^2]

y2 = @(x) rho0*cd0*(v_max^2)./(2*9.81*x)+2*9.81*K*x./(rho*sigma*(v_max^2))

plot(x2,y2(x2),'LineWidth',1.5)

y2(W_S)

% Passo 3 -> Take Off distance %

% la ground roll distance considerata � di 888 [m] ovvero circa 3000 [ft],
% ma la over 50 ft distance � maggiore e considerata come 888*1.7 [m] a cui
% corrisponde un valore di TOP pari a circa 230 [lbs/(ft^2)]
% e inoltre si considera un valore di clto pari a cl_max/1.21 come
% suggerito dal Roskam a pag 107

TOP = 230*0.45/(0.3048^2);           % [kg/m^2]
cl_to = cl_max/1.21;                

y3 = @(x) (x)./(TOP*cl_to);

plot(x2,y3(x2),'LineWidth',1.5)
axis([0 700 0 1])

Efficienza_max = 20;

% Passo 4 %
cd0_to = 0.075+0.02;
ROC = 3000*0.00508;     % [m/s]
y4 = @(x) ROC./sqrt(2*9.81*x./(rho0*sqrt(cd0_to/K)))+1/Efficienza_max;
plot(x2,y4(x2),'LineWidth',1.5)


% Passo 5 %
yline(1/Efficienza_max,'LineWidth',1.5)

% Punto di design
plot(W_S,y4(W_S),'*r','LineWidth',2)

legend('STALLO','CRUISE','TOP','ROC','CEILING','Design Point','Location','best')







