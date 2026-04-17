function dY = dinamica_senza_controllo(~, Y, p)
    % Estrazione stati
    S = max(Y(1), 0); % Cellule staminali ematopoietiche
    A = max(Y(2), 0); % Cellule progenitrici (sane)
    D = max(Y(3), 0); % Cellule differenziate (sane)
    L = max(Y(4), 0); % Cellule staminali leucemiche 
    T = max(Y(5), 0); % Cellule leucemiche differenziate (malate)
    
    % Variabili di popolazione per la competizione nello spazio midollare
    Z1 = S;
    Z2 = A + L;
    
    % Termine risposta immunitaria (Michaelis-Menten)
    % Se alpha = 0 (Modello Originale - Senza risposta immunitaria) il termine si annulla
    immune_response = (p.alpha * L) / (p.gamma + L);
    
    % Equazioni Differenziali (Eq. 2 del paper)
    dS = p.rho_S * S * (p.K1 - Z1) - p.delta_S * S;                   % Crescita logistica di S e differenziazione
    dA = p.delta_S * S + p.rho_A * A * (p.K2 - Z2) - p.delta_A * A;   % Differenziazione da S, crescita logistica di A (compete con L), differenziazione
    dD = p.delta_A * A - p.mu_D * D;                                  % Differenziazione da A, morte naturale/migrazione
    dL = p.rho_L * L * (p.K2 - Z2) - p.delta_L * L - immune_response; % Crescita logistica di L (compete con A), differenziazione, e risposta immunitaria
    dT = p.delta_L * L - p.mu_T * T;                                  % Differenziazione da L, morte naturale/migrazione
    
    dY = [dS; dA; dD; dL; dT];
end
