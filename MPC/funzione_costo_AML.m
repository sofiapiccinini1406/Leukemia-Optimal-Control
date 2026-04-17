function J = funzione_costo_AML(U_sequence, x_start, par, ind, Np, dt)
    % U_sequence: vettore dei controlli futuri proposti dall'ottimizzatore (lunghezza Np)
    % x_start: stato attuale del paziente da cui parte la predizione [S, A, D, L, T]
    
    a1 = ind(1); % Peso Farmaco
    a2 = ind(2); % Peso Leucemia
    
    J = 0;
    x_curr = x_start(:); % Assicura vettore colonna
    
    % Simuliamo in avanti per l'orizzonte predittivo Np
    for k = 1:Np
        u_k = U_sequence(k);
        
        % Integrazione Eulero (veloce per l'ottimizzazione)
        dx = dinamica_modello_interne(x_curr, u_k, par);
        x_next = x_curr + dx * dt;
        
        x_next = max(x_next, 0); % Vincolo non negatività
        L_val = x_next(4); % Estrazione Leucemia L
        
        % Calcolo Costo Istantaneo (Eq. 13 del paper: a1*u^2 + a2*L^2)
        step_cost = a1 * (u_k^2) + a2 * (L_val^2);
        
        J = J + step_cost * dt; % Integrale discreto
        x_curr = x_next; % Aggiorno stato
    end
end
