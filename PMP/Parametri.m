%% Parametri modello AML (Tabella 1 del Paper)
% Tassi di proliferazione
rho_S = 0.5;
rho_A = 0.43;
rho_L = 0.27;

% Tassi di differenziazione
delta_S = 0.14;
delta_A = 0.44;
delta_L = 0.05;

% Tassi di mortalità/migrazione
mu_D = 0.275;
mu_T = 0.3;

% Parametri risposta immunitaria
alpha = 0.015;
gamma = 0.1; 

% Capacità portanti
K1 = 1;
K2 = 1;

% Raggruppamento parametri in un vettore per passarli alle funzioni
par = [rho_S, rho_A, rho_L, delta_S, delta_A, delta_L, mu_D, mu_T, alpha, gamma, K1, K2];

%% Pesi per la funzione costo (Eq. 13 del Paper)
% J = integral( a1*u^2 + a2*L^2 )
a1 = 1;     % Peso costo chemioterapia
a2 = 1;     % Peso danno leucemia 

ind = [a1, a2]; 

%% Condizioni Iniziali (Coesistenza - Fig 4a)
S0 = 0.72;
A0 = 0.32;
D0 = 0.52;
L0 = 0.37;
T0 = 0.06; % Approssimati allo steady state

X_0 = [S0, A0, D0, L0, T0];

%% Parametri temporali
T = 50;                     % Durata simulazione
N = 1000;                   % Numero punti griglia
time = linspace(0, T, N);   
dt = time(2) - time(1);
