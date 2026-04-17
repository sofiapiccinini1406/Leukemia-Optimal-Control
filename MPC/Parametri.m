%% PARAMETRI MODELLO AML (Paper Tabella 1)
% Tassi di proliferazione
rho_S = 0.5;
rho_A = 0.43; 
rho_L = 0.27;

% Tassi di differenziazione
delta_S = 0.14; 
delta_A = 0.44; 
delta_L = 0.05;

% Tassi di mortalità
mu_D = 0.275; 
mu_T = 0.3;

% Immunità
alpha = 0.015;
gamma_val = 0.1;

% Capacità
K1 = 1; K2 = 1;

par = [rho_S, rho_A, rho_L, delta_S, delta_A, delta_L, mu_D, mu_T, alpha, gamma_val, K1, K2];

%% PESI FUNZIONE COSTO (J = a1*u^2 + a2*L^2)
% Usiamo il caso standard del paper (Fig 5b)
a1 = 0.1;  % Costo farmaco
a2 = 2;  % Costo leucemia
ind = [a1, a2];

%% PARAMETRI TEMPORALI E MPC
T_end = 50;                 % Durata totale simulazione
dt = 0.1;                   % Passo di discretizzazione (ogni 2.4 ore)
time = 0:dt:T_end;          % Vettore tempo totale
Nt = length(time);          % Numero passi totali

Np = 5;                     % Orizzonte di predizione
                            % Guarda avanti di 0.5 unità di tempo

% Limiti sul controllo u (Chemioterapia)
u_min = 0;                  % Minimo farmaco
u_max = 2.0;                % Massimo farmaco

%% CONDIZIONI INIZIALI (Coesistenza)
S0 = 0.72; A0 = 0.32; D0 = 0.52; L0 = 0.37; T0 = 0.06;
X0 = [S0, A0, D0, L0, T0];