%% MAIN MPC LEUCEMIA (Closed Loop)
clear; close all; clc;

%% Caricamento Parametri
Parametri; 
% Opzioni per l'ottimizzatore (comuni a entrambe le parti)
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'MaxIterations', 100);

%% PARTE 1: simulazione MPC - condizioni ideali
fprintf('--- PARTE 1: Simulazione MPC Ideale ---\n');

% Inizializzazione Simulazione Closed Loop
Np_ideal = Np; 

X_mpc_id = zeros(Nt, 5);       % Storia degli stati
U_mpc_id = zeros(Nt, 1);       % Storia dei controlli applicati
X_mpc_id(1,:) = X0;            % Stato iniziale

% Guess iniziale per l'ottimizzatore (vettore di zeri lungo Np)
u_guess_id = zeros(Np_ideal, 1);

% Ciclo MPC
for k = 1:Nt-1
    x_k = X_mpc_id(k, :)'; % Stato attuale del paziente
    
    % Definizione funzione costo per l'orizzonte corrente
    % Passiamo lo stato attuale x_k alla funzione costo
    cost_fun_id = @(U_seq) funzione_costo_AML(U_seq, x_k, par, ind, Np_ideal, dt);
    
    % Vincoli sui controlli per l'orizzonte (0 <= u <= u_max)
    lb_id = u_min * ones(Np_ideal, 1);
    ub_id = u_max * ones(Np_ideal, 1);
    
    % OTTIMIZZAZIONE (Trova la strategia migliore per i prossimi Np passi)
    [u_opt_seq, ~] = fmincon(cost_fun_id, u_guess_id, [],[],[],[], lb_id, ub_id, [], options);
    
    % Applichiamo SOLO il primo controllo calcolato
    u_applied = u_opt_seq(1);
    U_mpc_id(k) = u_applied;
    
    % Simuliamo il sistema "Reale" per un passo con questo controllo (uso ode45 per la massima precisione nella simulazione del paziente)
    t_span = [time(k), time(k+1)];
    [~, Y_step] = ode45(@(t,x) dinamica_modello_mpc(t, x, par, u_applied, dt), t_span, x_k);
    
    % Aggiorniamo lo stato "Reale" al passo successivo
    X_mpc_id(k+1, :) = Y_step(end, :);
    
    % Uso la soluzione vecchia come guess per il prossimo passo
    % Shiftiamo: [u2, u3, ..., uN, uN]
    u_guess_id = [u_opt_seq(2:end); u_opt_seq(end)];
    
    if mod(k, 50) == 0
        fprintf('Passo %d/%d completato (Leucemia: %.4f)\n', k, Nt, X_mpc_id(k+1,4));
    end
end

%% PARTE 2: simulazione CON DISTURBO - confronto open loop vs MPC
fprintf('--- PARTE 2: Stress Test (Ricaduta a t=15) ---\n');
% Aumentiamo l'orizzonte predittivo: Prima era 5. 

Np_dist = 20; % Portandolo a 20 (2 unità di tempo), l'MPC diventa meno "miope" e capisce che deve dare una dose alta all'inizio, simile all'Open Loop.

% Calcolo Piano OPEN LOOP (Strategia fissa calcolata a t=0)
fprintf('Calcolo strategia Open Loop (terapia fissa)...\n');
cost_fun_ol = @(U_seq) funzione_costo_AML(U_seq, X0, par, ind, Nt-1, dt);

% Guess iniziale un po' più intelligente
u_guess_ol = 0.1 * ones(Nt-1, 1);
lb_ol = u_min * ones(Nt-1, 1); 
ub_ol = u_max * ones(Nt-1, 1);

[U_open_loop, ~] = fmincon(cost_fun_ol, u_guess_ol, [],[],[],[], lb_ol, ub_ol, [], options);

% Simulazione pazienti
X_real_OL = zeros(Nt, 5); X_real_OL(1,:) = X0; % Open Loop
X_real_MPC = zeros(Nt, 5); X_real_MPC(1,:) = X0; % MPC
U_mpc_history = zeros(Nt-1, 1);

% Inizializzazione guess MPC (lungo Np)
u_guess_mpc = 0.1 * ones(Np_dist, 1);

for k = 1:Nt-1
    t_curr = time(k);
    
    u_ol_curr = U_open_loop(k);   % Strategia OPEN LOOP (Cieca)   
    x_k_mpc = X_real_MPC(k, :)';  % Strategia CLOSED LOOP (Reattiva)
    
    % Funzione costo MPC sull'orizzonte Np
    cost_mpc = @(U_seq) funzione_costo_AML(U_seq, x_k_mpc, par, ind, Np_dist, dt);
    lb_mpc = u_min * ones(Np_dist, 1); 
    ub_mpc = u_max * ones(Np_dist, 1);
    
    % Ottimizzazione rapida
    [u_seq_opt, ~] = fmincon(cost_mpc, u_guess_mpc, [],[],[],[], lb_mpc, ub_mpc, [], options);
    
    % Applica solo il primo controllo
    u_mpc_curr = u_seq_opt(1);
    U_mpc_history(k) = u_mpc_curr;
    
    u_guess_mpc = [u_seq_opt(2:end); u_seq_opt(end)]; % Shiftiamo la soluzione per il prossimo passo
    t_span = [time(k), time(k+1)]; % Evoluzione Reale
    
    % Simulazione Open Loop
    [~, Y_OL] = ode45(@(t,x) dinamica_modello_interne(x, u_ol_curr, par), t_span, X_real_OL(k,:));
    X_real_OL(k+1,:) = Y_OL(end,:);
    
    % Simulazione MPC
    [~, Y_MPC] = ode45(@(t,x) dinamica_modello_interne(x, u_mpc_curr, par), t_span, X_real_MPC(k,:));
    X_real_MPC(k+1,:) = Y_MPC(end,:);
    
    % Disturbo: al tempo t=15, simuliamo una ricaduta improvvisa (aumento L)
    if abs(t_curr - 15) < dt/2
        fprintf('-> DISTURBO: Ricaduta Improvvisa!\n');
        disturbo = 0.3;                                    % Aumento più marcato per rendere evidente il grafico
        X_real_OL(k+1, 4) = X_real_OL(k+1, 4) + disturbo;
        X_real_MPC(k+1, 4) = X_real_MPC(k+1, 4) + disturbo;
    end
    
    if mod(k, 50) == 0, fprintf('Step %d/%d...\n', k, Nt); end
end

%% Risultati Grafici
col_S = [0.85, 0.33, 0.10]; % Rosso
col_A = [1.00, 0.60, 0.20]; % Arancione
col_D = [0.92, 0.90, 0.15]; % Giallo
col_L = [0.05, 0.05, 0.40]; % Blu scuro
col_T = [0.39, 0.58, 0.93]; % Azzurro

% Figura 1: PLOT Dinamica delle Popolazioni
figure('Name', 'Dinamica MPC Closed Loop', 'NumberTitle', 'off');
hold on; grid on;

% Plot delle curve
p1 = plot(time, X_mpc_id(:,1), 'Color', col_S, 'LineWidth', 1.5); % S
p2 = plot(time, X_mpc_id(:,2), 'Color', col_A, 'LineWidth', 1.5); % A
p3 = plot(time, X_mpc_id(:,3), 'Color', col_D, 'LineWidth', 1.5); % D
p4 = plot(time, X_mpc_id(:,4), 'Color', col_L, 'LineWidth', 1.5); % L
p5 = plot(time, X_mpc_id(:,5), 'Color', col_T, 'LineWidth', 1.5); % T

ylabel('Popolazione Cellulare', 'FontSize', 12);
ylim([0 1]); % Scala fissa per le cellule
set(gca, 'YColor', 'k');

% ASSE DESTRO: Controllo
yyaxis right
% Plot del controllo (Tratteggiato Nero)
p6 = plot(time, U_mpc_id, 'k--', 'LineWidth', 2); 

ylabel('Dose Chemioterapia u(t)', 'FontSize', 12);
% ylim([0 u_max]); 
set(gca, 'YColor', 'k');
title('Dinamica Closed Loop (MPC)', 'FontSize', 14); xlabel('Tempo [giorni]', 'FontSize', 12); xlim([0 T_end]);
legend([p1, p2, p3, p4, p5, p6], {'S', 'A', 'D', 'L', 'T', 'u'}, 'Location', 'eastoutside', 'FontSize', 11);

% Figura 2: Stress Test (Disturbo)
figure('Name', 'Stress Test: Open Loop vs MPC', 'NumberTitle', 'off');

subplot(2,1,1); % Plot Superiore: Leucemia
% Area di "fallimento" Open Loop
x_fill = [time(time>=15), fliplr(time(time>=15))];
y1 = X_real_OL(time>=15, 4)';
y2 = X_real_MPC(time>=15, 4)';
fill(x_fill, [y1, fliplr(y2)], [1 0.9 0.9], 'EdgeColor', 'none'); hold on;

plot(time, X_real_OL(:,4), 'k--', 'LineWidth', 1.5); 
plot(time, X_real_MPC(:,4), 'Color', [0.05, 0.05, 0.40], 'LineWidth', 1.5);

ylabel('Popolazione Leucemica (L)', 'FontSize', 12); 
title('Effetto della Terapia Adattiva (MPC) vs Terapia Fissa', 'FontSize', 14);
legend('Area di Rischio (Open Loop)', 'Leucemia (Open Loop)', 'Leucemia (MPC)', 'Location', 'northeast');
grid on; ylim([0 1]); box on;

subplot(2,1,2); % Plot Inferiore: Dosaggio Farmaco 
plot(time(1:end-1), U_open_loop, 'm--', 'LineWidth', 1.5); hold on;
plot(time(1:end-1), U_mpc_history, 'Color', [0.05, 0.05, 0.4], 'LineWidth', 1.5);
xline(15, 'k--', 'LineWidth', 1.5);

ylabel('Dose Chemioterapia (u)', 'FontSize', 12); 
xlabel('Tempo [giorni]', 'FontSize', 12);
title('Adattamento del Dosaggio', 'FontSize', 14);
legend('Piano Originale (Pre-calcolato a t=0)', 'Reazione MPC (Real-time)');
grid on; box on;