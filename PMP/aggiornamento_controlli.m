function U_new = aggiornamento_controlli(X, P, U_old, N, par, ind, tipo_costo)
    U_new = zeros(N, 1);
    
    % Pesi costo
    a1 = ind(1);
    
    % Parametro di convergenza (omega)
    if strcmp(tipo_costo, 'bang-bang')
        omega = 0.9; % Il paper suggerisce valori alti (0.9) per il bang-bang per stabilità
        u_max = 0.5; % Limite superiore per Fig. 7
    else
        omega = 0.7; % Valore standard per il continuo
    end
    
    for k = 1:N

        L = X(k, 4);       % Stato Leucemia
        lambda4 = P(k, 4); % Co-stato Leucemia
        
        % Selezione del tipo di controllo
        if strcmp(tipo_costo, 'bang-bang')
            % CASO BANG-BANG (Costo Lineare)
            % Calcolo della Funzione di Commutazione
            phi = a1 - lambda4 * L;
            
            % Legge di controllo "On-Off"
            if phi < 0
                u_star = u_max; % Conviene curare (beneficio > costo)
            else
                u_star = 0;     % Non conviene curare
            end
            
        else
            % CASO CONTINUO (Costo Quadratico - DEFAULT)
            % Formula derivata da dH/du = 2*a1*u - lambda4*L = 0
            u_star = (lambda4 * L) / (2 * a1);
            
            % Vincolo di non negatività
            u_star = max(0, u_star); 
        end

        % Aggiornamento ponderato (Eq. 23 del Paper)
        if k > 1
             u_prev_iter = U_old(k);
             U_new(k) = omega * u_prev_iter + (1 - omega) * u_star;
        else
             U_new(k) = u_star;
        end
    end
end