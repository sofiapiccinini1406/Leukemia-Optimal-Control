function dx = dinamica_modello_interne(x, u, par)
    % Versione compatta della dinamica
    S=x(1); A=x(2); D=x(3); L=x(4); T=x(5);

    rho_S=par(1); rho_A=par(2); rho_L=par(3);
    delta_S=par(4); delta_A=par(5); delta_L=par(6);
    mu_D=par(7); mu_T=par(8);
    alpha=par(9); gamma_const=par(10);
    K1=par(11); K2=par(12);
    
    Z1 = S; Z2 = A + L;
    immune = (alpha * L) / (gamma_const + L);
    
    dx = zeros(5,1);
    dx(1) = rho_S * S * (K1 - Z1) - delta_S * S;
    dx(2) = delta_S * S + rho_A * A * (K2 - Z2) - delta_A * A;
    dx(3) = delta_A * A - mu_D * D;
    dx(4) = rho_L * L * (K2 - Z2) - delta_L * L - immune - u * L; % Controllo 
    dx(5) = delta_L * L - mu_T * T;
end