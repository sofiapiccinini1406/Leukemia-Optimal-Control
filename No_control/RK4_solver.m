%% RK4 integration on Acute Myeloid Leukaemia (AML) model
% Questo script simula il modello di AML modificato con risposta immunitaria (Eq. 2 del paper)

function [t, Y] = RK4_solver(odefun, tspan, y0, params)
    % Implementazione Runge-Kutta 4 standard
    h = tspan(2) - tspan(1);
    N = length(tspan);
    neq = length(y0);
    
    Y = zeros(N, neq);
    Y(1, :) = y0;
    
    for i = 1:N-1
        ti = tspan(i);
        yi = Y(i, :)';
        
        k1 = odefun(ti, yi, params);
        k2 = odefun(ti + 0.5*h, yi + 0.5*h*k1, params);
        k3 = odefun(ti + 0.5*h, yi + 0.5*h*k2, params);
        k4 = odefun(ti + h, yi + h*k3, params);
        
        ynew = yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        
        % Vincolo non-negatività semplice
        Y(i+1, :) = max(ynew, 0)';
    end
    t = tspan;
end
