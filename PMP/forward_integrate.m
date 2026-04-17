% Integra in avanti nel tempo le equazioni di stato del modello dinamico utilizzando il metodo di Runge–Kutta del 4° ordine (RK4).
% Dato lo stato iniziale e la sequenza dei controlli U, calcola la traiettoria completa
% degli stati X dal tempo t=0 al tempo t=T (ogni riga di X corrisponde a un tempo t)

function X = forward_integrate(x0, U, dt, N, par)
    X = zeros(N, 5); % 5 stati
    X(1,:) = x0; 

    for k = 1:N-1        
        u_k = U(k,:)';                                                     
        x_k = X(k,:)';
        
        k1 = dinamica_modello(x_k, u_k, par);
        k2 = dinamica_modello(x_k + 0.5*dt * k1, u_k, par);
        k3 = dinamica_modello(x_k + 0.5*dt * k2, u_k, par);
        k4 = dinamica_modello(x_k + dt * k3, u_k, par);
        
        X(k+1,:) = X(k,:) + (dt/6)*(k1' + 2*k2' + 2*k3' + k4');
        
        % Vincolo non negatività
        X(k+1,:) = max(X(k+1,:), 0);
    end
end
