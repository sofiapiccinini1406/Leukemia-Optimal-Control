% Questa funzione implementa il metodo di forward–backward sweep per risolvere un problema di controllo ottimo tramite il PMP
% L’algoritmo procede iterativamente:
%   1) integrazione in avanti delle equazioni di stato con i controlli attuali;
%   2) integrazione all’indietro delle equazioni dei co-stati;
%   3) aggiornamento dei controlli tramite le condizioni di ottimalità del PMP
%
% Le iterazioni continuano finché i controlli convergono (o viene raggiunto il numero massimo di iterazioni).
% La funzione restituisce le traiettorie ottimali di stati, co-stati e controlli.

function [X, P, U] = pmp_forward_backward_solver(X_0, N, max_iter, par, ind, dt, tipo_costo)
    
% Se non specificato, assumiamo costo quadratico (default)
    if nargin < 7
        tipo_costo = 'quadratico';
    end

    % Inizializzazione controllo (Guess iniziale u=0 come nel paper)
    U = zeros(N, 1);    
    
    % Tolleranza per convergenza
    tol = 1e-3;
    
    for iter = 1:max_iter
        U_old = U;
        
        % 1) Integrazione in avanti (Stati)
        X = forward_integrate(X_0, U, dt, N, par);
        
        % 2) Integrazione all'indietro (Co-stati)
        P = backward_integrate(X, U, dt, N, par, ind, tipo_costo);
        
        % 3) Aggiornamento controlli (PMP + media ponderata)
        U = aggiornamento_controlli(X, P, U, N, par, ind, tipo_costo);
        
        % Verifica convergenza (Eq. 21 del paper)
        diff_u = sum(abs(U - U_old)) / (sum(abs(U)) + 1e-10);
        
        if diff_u < tol
            fprintf('Convergenza raggiunta all''iterazione %d. Errore: %e\n', iter, diff_u);
            break;
        end
        if iter == max_iter
            fprintf('Attenzione: Max iterazioni raggiunte senza convergenza completa.\n');
        end
    end
end
