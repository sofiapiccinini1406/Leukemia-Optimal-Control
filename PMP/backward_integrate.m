% backward_integrate
% Integra all’indietro nel tempo le equazioni dei co-stati derivate dal Principio del Massimo di Pontryagin (PMP).
% Utilizza il metodo di Runge–Kutta del 4° ordine (RK4) e impone le condizioni finali dei 
% co-stati a t = T, restituendo l’intera traiettoria dei co-stati P(t).
function P = backward_integrate(X, U, dt, N, par, ind, tipo_costo)
    P = zeros(N, 5);        % 5 co-stati
    P(end,:) = [0 0 0 0 0]; % Condizione di trasversalità (lambda(tf) = 0)
    
    for k = N:-1:2                                                                                                                                    
        x_k = X(k,:)';                                                      
        p_k = P(k,:)';                                                                                                                         
        u_k = U(k,:)';                                                      
        
        k1 = co_stati(x_k, p_k, u_k, par, ind, tipo_costo);
        k2 = co_stati(x_k, p_k - 0.5*dt*k1, u_k, par, ind, tipo_costo);
        k3 = co_stati(x_k, p_k - 0.5*dt*k2, u_k, par, ind, tipo_costo);
        k4 = co_stati(x_k, p_k - dt*k3, u_k, par, ind, tipo_costo);
        
        % Integrazione all'indietro: P(k-1) = P(k) - dt*media_pendenze
        P(k-1,:) = P(k,:) - (dt/6)*(k1' + 2*k2' + 2*k3' + k4');
    end
end