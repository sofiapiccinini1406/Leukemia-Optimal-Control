clear; close all; clc;

%% 1. DEFINIZIONE PARAMETRI DA TABELLA 1
rho_S   = 0.5;   % Differenziazione di S in A
rho_A   = 0.43;  % Differenziazione di A in D
rho_L   = 0.27;  % Differenziazione di L in T
delta_S = 0.14;  % Proliferazione di S (Staminali sane)
delta_A = 0.44;  % Proliferazione di A (Progenitori sani)
delta_L = 0.05;  % Proliferazione di L (Staminali leucemiche)
mu_D    = 0.275; % Migrazione di D nel flusso sanguigno
mu_T    = 0.3;   % Migrazione di T nel flusso sanguigno
K1      = 1;     % Capacità portante compartimento con S
K2      = 1;     % Capacità portante compartimento con A e L

% Parametri Risposta Immunitaria (Tabella 1)
alpha_tab = 0.015; % Risposta immunitaria
gamma_tab = 0.01;  % Saturazione risposta immunitaria

% Raggruppo i parametri fissi in una struct
par_base.rho_S = rho_S;      par_base.rho_A = rho_A;        par_base.rho_L = rho_L;      
par_base.delta_S = delta_S;  par_base.delta_A = delta_A;    par_base.delta_L = delta_L;
par_base.mu_D = mu_D;        par_base.mu_T = mu_T;
par_base.K1 = K1;            par_base.K2 = K2;

%% 2. PARTE (A): MODELLO ORIGINALE
% Riferimento: Figura 2(a, b, c, d) del paper
% Il modello originale NON ha risposta immunitaria. Lo otteniamo ponendo alpha = 0

fprintf('Simulazione (a): Modello Originale \n');

% Parametri specifici per il grafico (a)
par_a = par_base;
par_a.alpha = 0;          % Nessuna risposta immunitaria
par_a.gamma = gamma_tab;  % Irrilevante con alpha=0

% Condizioni Iniziali (didascalia Fig. 2a: Coesistenza)
% [S, A, D, L, T]
X0_a = [0.1, 0, 0, 0.1, 0];  % coesistenza
X0_b = [0.5, 0, 0, 1e-3, 0]; % S parte alto, L parte molto basso
X0_c = [0.1, 0, 0, 0, 0];    % stato stazionario sano
X0_d = [0, 0, 0, 0.1, 0];    % stato stazionario malato

% Tempo
tspan = 0:0.1:50; % arriva a 50
tspan_b = 0:0.1:150; % arriva a 150

% Risoluzione
[T_a, Y_a] = RK4_solver(@dinamica_senza_controllo, tspan, X0_a, par_a);
[T_b, Y_b] = RK4_solver(@dinamica_senza_controllo, tspan_b, X0_b, par_a);
[T_c, Y_c] = RK4_solver(@dinamica_senza_controllo, tspan, X0_c, par_a);
[T_d, Y_d] = RK4_solver(@dinamica_senza_controllo, tspan, X0_d, par_a);

% FIGURA 1: Grafico Modello originale (a)
figure('Name', 'Modello (a): Originale', 'NumberTitle', 'off');
plot_aml_results(T_a, Y_a, 'Nessuna Risposta Immunitaria - Stato Stazionario di Coesistenza');

% FIGURA 2: Grafico Modello originale (b)
figure('Name', 'Modello (b): popolazione leucemica iniziale bassa', 'NumberTitle', 'off');
plot_aml_results(T_b, Y_b, 'Nessuna Risposta Immunitaria - Stato Stazionario di Coesistenza');

% FIGURA 3: Griglia Comparativa 2x2
figure('Name', 'Confronto Evoluzioni Modello Originale', 'NumberTitle', 'off');

% In alto a sinistra: Primo grafico
subplot(2, 2, 1); plot_aml_results(T_a, Y_a, '(a) Coesistenza');

% In alto a destra: Secondo grafico
% non è possibile raggiungere un stato stazionario sano non leucemico in presenza anche di piccole popolazioni di cellule staminali leucemiche
subplot(2, 2, 2); plot_aml_results(T_b, Y_b, '(b) Coesistenza (popolazione leucemica iniziale bassa)');

% In basso a sinistra: Terzo grafico
subplot(2, 2, 3); plot_aml_results(T_c, Y_c, '(c) Stato stazionario Sano');

% In basso a destra: Quarto grafico
subplot(2, 2, 4); plot_aml_results(T_d, Y_d, '(d) Stato stazionario Leucemico');

%% 3. PARTE (B): MODELLO MODIFICATO
% Riferimento: Figura 4(a, b, c, d) del paper.
% Il modello include la risposta immunitaria (efficace per L piccolo e inefficace per L grande)
% NOTA: La didascalia della Fig. 4 specifica gamma = 0.1 (diverso dalla Tabella 1)

fprintf('Simulazione (b): Modello Modificato con Risposta Immunitaria \n');

% Parametri specifici per il grafico (b)
par_b = par_base;
par_b.alpha = alpha_tab;  % 0.015 (da Tabella 1)
par_b.gamma = 0.1;        % Valore specifico per Fig. 4

% Condizioni Iniziali uguali alle precedenti X0_a, X0_b, X0_c, X0_d 

% Risoluzione
[T_a, Y_a] = RK4_solver(@dinamica_senza_controllo, tspan_b, X0_a, par_b);
[T_b, Y_b] = RK4_solver(@dinamica_senza_controllo, tspan, X0_b, par_b);
[T_c, Y_c] = RK4_solver(@dinamica_senza_controllo, tspan, X0_c, par_b);
[T_d, Y_d] = RK4_solver(@dinamica_senza_controllo, tspan, X0_d, par_b);

% FIGURA 4: Grafico Modello modificato (a)
figure('Name', 'Modello (a): Modificato', 'NumberTitle', 'off');
plot_aml_results(T_a, Y_a, 'Modello Modificato con Risposta Immunitaria');

% FIGURA 5: Grafico Modello modificato (b)
figure('Name', 'Modello (b): Modificato', 'NumberTitle', 'off');
plot_aml_results(T_b, Y_b, 'Modello Modificato con Risposta Immunitaria');

% FIGURA 6: Griglia Comparativa 2x2
figure('Name', 'Confronto Evoluzioni Modello Modificato', 'NumberTitle', 'off');

% In alto a sinistra: Primo grafico
% impiega più tempo per approcciare lo stato stazionario rispetto al modello originale 
subplot(2, 2, 1); plot_aml_results(T_a, Y_a, '(a) Coesistenza');

% In alto a destra: Secondo grafico
% con la risposta immunitaria una piccola popolazione di cellule staminali leucemiche non sopravvive
subplot(2, 2, 2); plot_aml_results(T_b, Y_b, '(b) Popolazione leucemica iniziale bassa - non sopravvive');

% In basso a sinistra: Terzo grafico
subplot(2, 2, 3); plot_aml_results(T_c, Y_c, '(c) Stato stazionario Sano');

% In basso a destra: Quarto grafico
subplot(2, 2, 4); plot_aml_results(T_d, Y_d, '(d) Stato stazionario Leucemico');

%  FUNZIONE LOCALE
function plot_aml_results(t, Y, title_str)
    hold on; grid on;
    % Colori e stili
    col_S = [0.85, 0.33, 0.10]; % S: Rosso
    col_A = [1.00, 0.60, 0.20]; % A: Arancione
    col_D = [0.92, 0.90, 0.15]; % D: Giallo
    col_L = [0.05, 0.05, 0.40]; % L: Blu scuro   
    col_T = [0.39, 0.58, 0.93]; % T: Azzurro chiaro

    plot(t, Y(:,1), '-', 'Color', col_S, 'LineWidth', 1.5, 'DisplayName', 'S');
    plot(t, Y(:,2), '-', 'Color', col_A, 'LineWidth', 1.5, 'DisplayName', 'A');
    plot(t, Y(:,3), '-', 'Color', col_D, 'LineWidth', 1.5, 'DisplayName', 'D');
    plot(t, Y(:,4), '-', 'Color', col_L, 'LineWidth', 1.5, 'DisplayName', 'L');
    plot(t, Y(:,5), '-', 'Color', col_T, 'LineWidth', 1.5, 'DisplayName', 'T');
    
    xlabel('Tempo [giorni]');
    ylabel('Popolazione Cellulare');
    title(title_str);
    legend('Location', 'best');
    ylim([0 1]);
    box on;
end