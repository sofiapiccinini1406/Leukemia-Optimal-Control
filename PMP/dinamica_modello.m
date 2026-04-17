function dx = dinamica_modello(x, u, par)
    % Stati
    S = max(x(1),0); A = max(x(2),0); D = max(x(3),0); 
    L = max(x(4),0); T = max(x(5),0);
    
    % Parametri
    rho_S=par(1); rho_A=par(2); rho_L=par(3);
    delta_S=par(4); delta_A=par(5); delta_L=par(6);
    mu_D=par(7); mu_T=par(8);
    alpha=par(9); gamma_const=par(10);
    K1=par(11); K2=par(12);
    
    % Controllo
    chemo = u(1); 

    % Variabili accoppiate
    Z1 = S;
    Z2 = A + L;
    
    % Risposta immunitaria
    immune = (alpha * L) / (gamma_const + L);
    
    dx = zeros(5,1);
    
    dx(1) = rho_S * S * (K1 - Z1) - delta_S * S;                      % dS/dt
    dx(2) = delta_S * S + rho_A * A * (K2 - Z2) - delta_A * A;        % dA/dt
    dx(3) = delta_A * A - mu_D * D;                                   % dD/dt 
    dx(4) = rho_L * L * (K2 - Z2) - delta_L * L - immune - chemo * L; % dL/dt (con immunità e controllo -u*L)
    dx(5) = delta_L * L - mu_T * T;                                   % dT/dt
end

