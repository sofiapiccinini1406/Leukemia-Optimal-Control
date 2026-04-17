function dp = co_stati(x, p, u, par, ind, tipo_costo)
    % Stati
    S = x(1); A = x(2); D = x(3); L = x(4); T = x(5);
    
    % Co-stati (Lambda)
    l1 = p(1); l2 = p(2); l3 = p(3); l4 = p(4); l5 = p(5);
    
    % Controllo
    chemo = u(1);
    
    % Parametri e Pesi
    rho_S=par(1); rho_A=par(2); rho_L=par(3);
    delta_S=par(4); delta_A=par(5); delta_L=par(6);
    mu_D=par(7); mu_T=par(8);
    alpha=par(9); gamma_const=par(10);
    
    a1 = ind(1); % Peso controllo
    a2 = ind(2); % Peso leucemia
    
    % Eq. 15 del Paper: dLambda/dt    
    dp = zeros(5,1);
    
    dp(1) = 2*S*l1*rho_S + delta_S*l1 - delta_S*l2 - l1*rho_S;                           % dLambda1/dt
    dp(2) = 2*A*l2*rho_A + L*l2*rho_A + L*l4*rho_L + delta_A*l2 - delta_A*l3 - l2*rho_A; % dLambda2/dt
    dp(3) = mu_D * l3;                                                                   % dLambda3/dt  
    term_immune = (alpha * gamma_const * l4) / ((gamma_const + L)^2);                    % dLambda4/dt

    if strcmp(tipo_costo, 'quadratico')
        % CASO CONTINUO (Fig 5/6)
        dL_term = -2*a2*L; 
    else 
        % CASO BANG-BANG (Fig 7)
        dL_term = -a2;
    end
    dp(4) = dL_term + rho_A*A*l2 + l4*rho_L*A + 2*rho_L*L*l4 - l4*rho_L + l4*delta_L + term_immune + l4*chemo - delta_L*l5;
    dp(5) = mu_T * l5;                                                                    % dLambda5/dt
end

