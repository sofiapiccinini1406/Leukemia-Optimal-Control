%% MAIN PMP CONTROLLO CONTINUO: Analisi Controllo Ottimo Leucemia (Figure 5 e 6 paper)
clear; close all; clc;

%% INIZIALIZZAZIONE PARAMETRI
Parametri;

T = 50; 
N = 1000; 
time = linspace(0, T, N);
dt = time(2) - time(1);

%% PARTE 1: FIGURA 5 (Confronto Pesi a1, a2)
% Confronto scenari di controllo ottimo al variare dei pesi (a1, a2)
% SCENARIO A: stato stazionario di coesistenza senza controllo applicato
fprintf('Calcolo Scenario A (Nessun Controllo)...\n');
X_nc = zeros(N, 5); X_nc(1,:) = X_0;
for k=1:N-1
    % Integrazione manuale dinamica con u=0
    k1 = dinamica_modello(X_nc(k,:)', 0, par);
    k2 = dinamica_modello(X_nc(k,:)' + 0.5*dt*k1, 0, par);
    k3 = dinamica_modello(X_nc(k,:)' + 0.5*dt*k2, 0, par);
    k4 = dinamica_modello(X_nc(k,:)' + dt*k3, 0, par);
    X_nc(k+1,:) = X_nc(k,:) + (dt/6)*(k1'+2*k2'+2*k3'+k4');
    X_nc(k+1,:) = max(X_nc(k+1,:), 0);
end
U_nc = zeros(N,1); 

% SCENARIO B: Standard - peso uguale (a1=1, a2=1)
fprintf('Calcolo Scenario B (Standard: a1=1, a2=1)...\n');
ind_b = [1, 1]; 
[X_b, ~, U_b] = pmp_forward_backward_solver(X_0, N, 500, par, ind_b, dt);
J_b = calcola_J_quadratico(U_b, X_b(:,4), ind_b, time);
fprintf('      -> J = %.4f\n', J_b);

% SCENARIO C: Controllo aggressivo (a1=0.1, a2=1)
fprintf('Calcolo Scenario C (Aggressivo: a1=0.1)...\n');
ind_c = [0.1, 1]; 
[X_c, ~, U_c] = pmp_forward_backward_solver(X_0, N, 500, par, ind_c, dt);
J_c = calcola_J_quadratico(U_c, X_c(:,4), ind_c, time);
fprintf('      -> J = %.4f\n', J_c);

% SCENARIO D: Controllo conservativo (a1=1, a2=0.1)
fprintf('Calcolo Scenario D (Conservativo: a2=0.1)...\n');
ind_d = [1, 0.1]; 
[X_d, ~, U_d] = pmp_forward_backward_solver(X_0, N, 500, par, ind_d, dt);
J_d = calcola_J_quadratico(U_d, X_d(:,4), ind_d, time);
fprintf('      -> J = %.4f\n', J_d);

% PLOT FIGURA 5 (Griglia 2x2)
figure('Name', 'Confronto Strategie Ottime', 'NumberTitle', 'off');

subplot(2,2,1); plot_scenario_custom(time, X_nc, U_nc, '(a) Stato Stazionario di Coesistenza Senza Controllo', true); 
subplot(2,2,2); plot_scenario_custom(time, X_b, U_b, '(b) Standard - peso uguale [a1=1, a2=1]', true);
subplot(2,2,3); plot_scenario_custom(time, X_c, U_c, '(c) Controllo aggressivo [a1=0.1, a2=1]', true);
subplot(2,2,4); plot_scenario_custom(time, X_d, U_d, '(d) Controllo conservativo [a1=1, a2=0.1]', true);

%% Confronto tra modello senza controllo e modello controllato - plot dettaglio scenario B 
% Analisi delle popolazioni per il caso standard
% Nuova Condizione Iniziale (Sano alto, Leucemia bassa)
X_0_dyn = [0.7, 0.05, 0.05, 0.07, 0]; 
% X_0_dyn = X_0;

% Calcolo Dinamica NO CONTROL (u=0)
fprintf('Ricalcolo Scenario No-Control per grafico dinamico...\n');
X_nc_dyn = zeros(N, 5); X_nc_dyn(1,:) = X_0_dyn;
for k=1:N-1
    k1 = dinamica_modello(X_nc_dyn(k,:)', 0, par); % u=0
    k2 = dinamica_modello(X_nc_dyn(k,:)' + 0.5*dt*k1, 0, par);
    k3 = dinamica_modello(X_nc_dyn(k,:)' + 0.5*dt*k2, 0, par);
    k4 = dinamica_modello(X_nc_dyn(k,:)' + dt*k3, 0, par);
    X_nc_dyn(k+1,:) = X_nc_dyn(k,:) + (dt/6)*(k1'+2*k2'+2*k3'+k4');
    X_nc_dyn(k+1,:) = max(X_nc_dyn(k+1,:), 0);
end

% Calcolo Dinamica OPTIMAL CONTROL (PMP)
fprintf('Ricalcolo Scenario Ottimo per grafico dinamico...\n');
ind_b = [1, 1]; % Pesi standard
% Rilancio il solver partendo dalla nuova condizione iniziale
[X_b_dyn, ~, U_b_dyn] = pmp_forward_backward_solver(X_0_dyn, N, 500, par, ind_b, dt);
figure('Name', 'Confronto Diretto: No Control vs Optimal Control', 'NumberTitle', 'off');

labels = ["Cellule staminali Sane (S)", "Cellule progenitrici Sane (A)", "Cellule differenziate Sane (D)", "Cellule staminali Leucemiche (L)", "Cellule Leucemiche differenziate (T)"];
% Definizione Colori
col_S = [0.85, 0.33, 0.10]; % Rosso
col_A = [1.00, 0.60, 0.20]; % Arancione
col_D = [0.92, 0.90, 0.15]; % Giallo
col_L = [0.05, 0.05, 0.40]; % Blu scuro 
col_T = [0.39, 0.58, 0.93]; % Azzurro
colors_cell = {col_S, col_A, col_D, col_L, col_T};

% Plot dei 5 stati
for i = 1:5
    subplot(2,3,i);
    % No Control - Grigio tratteggiato
    plot(time, X_nc_dyn(:,i), 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 2); hold on;
    plot(time, X_b_dyn(:,i), 'Color', colors_cell{i}, 'LineWidth', 2);
    title(labels(i)); xlabel('Tempo [giorni]'); ylabel('Popolazione Cellulare');
    legend('No Control', 'Controllo Continuo', 'Location', 'best');
    grid on; ylim([0 1]); xlim([0 50]); box on;
end

% Plot del Controllo u(t)
subplot(2,3,6);
plot(time, U_b, 'k--', 'LineWidth', 2);
title('Controllo Ottimo u(t)');
xlabel('Tempo [giorni]'); ylabel('Intensità');
xlim([0 50]);
grid on; box on;

%% PARTE 2: FIGURA 6 (Confronto Strategie a Parità di Dosaggio)
% Confronto tra controllo ottimo e strategie euristiche a parità di dosaggio
% Usiamo lo Scenario B come Benchmark (Ottimo)
U_opt = U_b;
X_opt = X_b;
J_opt = J_b;
ind_std = ind_b; % Pesi standard per il confronto [1, 1]

% Calcolo del "Budget" di farmaco (Dosaggio Totale usato dall'ottimo)
Dosaggio_Totale = trapz(time, U_opt);
fprintf('Dosaggio Totale Disponibile: %.4f\n', Dosaggio_Totale);

% Creazione Strategie Alternative (Euristiche)
% b) Strategia COSTANTE
val_costante = Dosaggio_Totale / T;
U_const = ones(N, 1) * val_costante;
X_const = forward_integrate(X_0, U_const, dt, N, par);
J_const = calcola_J_quadratico(U_const, X_const(:,4), ind_std, time);

% c) Strategia CICLO UNICO (Short Cycle)
T_cycle1 = 5; 
val_cycle1 = Dosaggio_Totale / T_cycle1;
U_cycle1 = zeros(N, 1);
U_cycle1(time <= T_cycle1) = val_cycle1;
X_cycle1 = forward_integrate(X_0, U_cycle1, dt, N, par);
J_cycle1 = calcola_J_quadratico(U_cycle1, X_cycle1(:,4), ind_std, time);

% d) Strategia DUE CICLI (Two Cycles)
% Lunghezza cicli: [0-20] e [30-50]
t_start1 = 0;  t_end1 = 20;
t_start2 = 30; t_end2 = 50;
durata_totale_on = (t_end1 - t_start1) + (t_end2 - t_start2);

% L'intensità deve essere più bassa perché spalmata su più tempo
val_cycle2 = Dosaggio_Totale / durata_totale_on; 

U_cycle2 = zeros(N, 1);
% Primo ciclo (0 -> 20)
U_cycle2(time >= t_start1 & time <= t_end1) = val_cycle2;
% Secondo ciclo (30 -> 50)
U_cycle2(time >= t_start2 & time <= t_end2) = val_cycle2; 

X_cycle2 = forward_integrate(X_0, U_cycle2, dt, N, par);
J_cycle2 = calcola_J_quadratico(U_cycle2, X_cycle2(:,4), ind_std, time);

% PLOT FIGURA 6 (Griglia 2x2)
figure('Name', 'Confronto Efficienza Strategie', 'NumberTitle', 'off');

subplot(2,2,1); plot_scenario_custom(time, X_opt, U_opt, sprintf('(a) Controllo Ottimo\nJ = %.4f', J_opt), true);
subplot(2,2,2); plot_scenario_custom(time, X_const, U_const, sprintf('(b) Controllo a Tasso costante\nJ = %.4f', J_const), true);
subplot(2,2,3); plot_scenario_custom(time, X_cycle1, U_cycle1, sprintf('(c) Controllo a Tasso più alto (Ciclo Breve)\nJ = %.4f', J_cycle1), true);
subplot(2,2,4); plot_scenario_custom(time, X_cycle2, U_cycle2, sprintf('(d) Controllo in due cicli\nJ = %.4f', J_cycle2), true);

%% PARTE 3: FIGURA 7 (Controllo Bang-Bang)
tipo_simulazione = 'bang-bang'; % Passiamo il flag 'bang-bang' al solver
u_max_bb = 0.5; % Valore usato nella Fig 7

% d) CASO 1: Farmaco "Economico" (a1=0.1, a2=1)
% Poiché costa poco (0.1), sarà un trattamento LUNGO (fino a t=8)
fprintf('Calcolo Caso Bang-Bang 1 (Farmaco Economico: a1=0.1, a2=1)...\n');
ind_bb1 = [0.1, 1];
[X_bb1, ~, U_bb1] = pmp_forward_backward_solver(X_0, N, 2000, par, ind_bb1, dt, tipo_simulazione);
J_bb1 = calcola_J_lineare(U_bb1, X_bb1(:,4), ind_bb1, time);

% e) CASO 2: Bilanciato (a1=1, a2=1) 
% Trattamento MEDIO (fino a t=4)
fprintf('Calcolo Caso Bang-Bang 2 (Bilanciato: a1=1, a2=1)...\n');
ind_bb2 = [1, 1];
[X_bb2, ~, U_bb2] = pmp_forward_backward_solver(X_0, N, 2000, par, ind_bb2, dt, tipo_simulazione);
J_bb2 = calcola_J_lineare(U_bb2, X_bb2(:,4), ind_bb2, time);

% f) CASO 3: Farmaco "Costoso" / Leucemia Irrilevante (a1=1, a2=0.1) 
% Poiché il danno della malattia pesa poco (0.1), NON conviene curare - trattamento NULLO 
fprintf('Calcolo Caso Bang-Bang 3 (Farmaco Costoso: a1=1, a2=0.1)...\n');
ind_bb3 = [1, 0.1];
[X_bb3, ~, U_bb3] = pmp_forward_backward_solver(X_0, N, 2000, par, ind_bb3, dt, tipo_simulazione);
J_bb3 = calcola_J_lineare(U_bb3, X_bb3(:,4), ind_bb3, time);

% PLOT FIGURA 7 (Griglia 1x3)
figure('Name', 'Controllo Ottimo Bang-Bang', 'NumberTitle', 'off');

subplot(1,3,1); plot_bang_bang(time, X_bb1, U_bb1, u_max_bb, sprintf('Farmaco Economico (a_1 = 0.1, a_2 = 1)\n(Trattamento Lungo)'));
subplot(1,3,2); plot_bang_bang(time, X_bb2, U_bb2, u_max_bb, sprintf('Bilanciato (a_1 = 1, a_2 = 1)\n(Trattamento Medio)'));
subplot(1,3,3); plot_bang_bang(time, X_bb3, U_bb3, u_max_bb, sprintf('Farmaco Costoso (a_1 = 1, a_2 = 0.1)\n(Trattamento Minimo)'));

%% Funzioni Ausiliarie Locali

function J = calcola_J_quadratico(U, L, ind, t)
    % Costo Quadratico (Fig 5 e 6): a1*u^2 + a2*L^2
    integrand = ind(1)*(U.^2) + ind(2)*(L.^2);
    J = trapz(t, integrand);
end

function J = calcola_J_lineare(U, L, ind, t)
    % Costo Lineare (Fig 7 - Bang Bang): a1*u + a2*L
    integrand = ind(1)*U + ind(2)*L;
    J = trapz(t, integrand);
end

function plot_scenario_custom(t, X, U, title_str, plot_control)
    hold on; grid on;
    % Colori richiesti
    col_S = [0.85, 0.33, 0.10];
    col_A = [1.00, 0.60, 0.20]; 
    col_D = [0.92, 0.90, 0.15];
    col_L = [0.05, 0.05, 0.40]; 
    col_T = [0.39, 0.58, 0.93]; 
    
    plot(t, X(:,1), 'Color', col_S, 'LineWidth', 1.5);
    plot(t, X(:,2), 'Color', col_A, 'LineWidth', 1.5);
    plot(t, X(:,3), 'Color', col_D, 'LineWidth', 1.5);
    plot(t, X(:,4), 'Color', col_L, 'LineWidth', 1.5);
    plot(t, X(:,5), 'Color', col_T, 'LineWidth', 1.5);
    
    if plot_control
        plot(t, U, 'k--', 'LineWidth', 1.5); % Controllo nero tratteggiato
    end
    
    title(title_str);
    if contains(title_str, 'Standard') || contains(title_str, 'Constant')
        if plot_control
             legend('S','A','D','L','T','u*','Location','best');
        else
             legend('S','A','D','L','T','Location','best');
        end
    end
    
    xlabel('Tempo [giorni]');
    ylabel('Popolazione Cellulare');
    ylim([0 1]); xlim([0 50]);
    box on;
end

function plot_bang_bang(t, X, U, u_max, title_str)
    hold on; grid on;
    
    % Definizione Colori
    col_S = [0.85, 0.33, 0.10]; % Rosso 
    col_A = [1.00, 0.60, 0.20]; % Arancione
    col_D = [0.92, 0.90, 0.15]; % Giallo
    col_L = [0.05, 0.05, 0.40]; % Blu scuro
    col_T = [0.39, 0.58, 0.93]; % Azzurro
    
    % Plot degli Stati
    p1 = plot(t, X(:,1), 'Color', col_S, 'LineWidth', 1.5);
    p2 = plot(t, X(:,2), 'Color', col_A, 'LineWidth', 1.5);
    p3 = plot(t, X(:,3), 'Color', col_D, 'LineWidth', 1.5);
    p4 = plot(t, X(:,4), 'Color', col_L, 'LineWidth', 1.5); 
    p5 = plot(t, X(:,5), 'Color', col_T, 'LineWidth', 1.5);
    
    % Plot del Controllo (Tratteggiato nero)
    p6 = plot(t, U, 'k--', 'LineWidth', 1.5);
    
    title(title_str, 'FontSize', 10);
    xlabel('Tempo [giorni]');
    ylabel('Popolazione cellulare / Controllo');
    ylim([0 1]);      % Gli stati vanno da 0 a 1 (normalizzati)
    xlim([0 50]);
    
    % Legenda completa
    legend([p1 p2 p3 p4 p5 p6], {'S', 'A', 'D', 'L', 'T', 'u*'}, 'Location', 'east', 'FontSize', 8);
    box on;
end

%% PARTE 4: Simulazione "Drug Holidays" (Protocollo Clinico Bang-Bang Multi-Ciclo)
fprintf('Calcolo Strategia Drug Holidays (Cicli ON/OFF)...\n');

U_holidays = zeros(N, 1);
u_max_bb = 0.5; % Dose massima come nel bang-bang

% Definiamo la durata del ciclo terapeutico
giorni_ON = 2;   % Giorni di chemioterapia al massimo dosaggio
giorni_OFF = 10; % Giorni di "vacanza" (sospensione totale)
periodo = giorni_ON + giorni_OFF;

% Creiamo il segnale Bang-Bang periodico
for i = 1:N
    t_corrente = time(i);
    giorno_nel_ciclo = mod(t_corrente, periodo);
    if giorno_nel_ciclo <= giorni_ON
        U_holidays(i) = u_max_bb; % Farmaco acceso al massimo
    else
        U_holidays(i) = 0;        % Sospensione (Drug Holiday)
    end
end

% Integriamo la dinamica con il nostro protocollo Drug Holidays
X_holidays = zeros(N, 5); X_holidays(1,:) = X_0;
for k=1:N-1
    k1 = dinamica_modello(X_holidays(k,:)', U_holidays(k), par);
    k2 = dinamica_modello(X_holidays(k,:)' + 0.5*dt*k1, U_holidays(k), par);
    k3 = dinamica_modello(X_holidays(k,:)' + 0.5*dt*k2, U_holidays(k), par);
    k4 = dinamica_modello(X_holidays(k,:)' + dt*k3, U_holidays(k), par);
    X_holidays(k+1,:) = X_holidays(k,:) + (dt/6)*(k1'+2*k2'+2*k3'+k4');
    X_holidays(k+1,:) = max(X_holidays(k+1,:), 0);
end

% Calcoliamo il costo di questa strategia 
J_holidays = calcola_J_lineare(U_holidays, X_holidays(:,4), [1, 1], time);
fprintf('      -> J Drug Holidays = %.4f\n', J_holidays);

% PLOT FIGURA DRUG HOLIDAYS
figure('Name', 'Simulazione Clinica: Drug Holidays', 'NumberTitle', 'off');
plot_bang_bang(time, X_holidays, U_holidays, u_max_bb, sprintf('Protocollo Drug Holidays (ON: %d gg, OFF: %d gg)', giorni_ON, giorni_OFF));