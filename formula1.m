function risk = formula1(psi1,psi2,lambda,mu1,mustar,F1,Fstar,tau)
%FORMULA Summary of this function goes here
%   Detailed explanation goes here
    zeta = mu1 / mustar;
    lambda_bar = lambda / mustar^2;
    xi = complex(0,sqrt(psi1*psi2*lambda_bar));
    syms nu1 nu2
    eq1 = nu1+psi1/((xi+nu2+(zeta^2*nu2)/(1-zeta^2*nu1*nu2)));
    eq2 = nu2+psi2/((xi+nu1+(zeta^2*nu1)/(1-zeta^2*nu1*nu2)));
    rad1 = sqrt(psi1 / (lambda_bar * psi2));
    rad2 = sqrt(psi2 / (lambda_bar * psi1));
    init_param1 = [complex(-rad1,0) complex(rad1,rad1)];
    init_param2 = [complex(-rad2,0) complex(rad2,rad2)];
    result = vpasolve([eq1,eq2],[nu1,nu2],[init_param1; init_param2]);
    nu1 = double(result.nu1);
    nu2 = double(result.nu2);
    chi = nu1 * nu2;
    %% Generalization error formula
    E0 = -chi^5*zeta^6+3*chi^4*zeta^4+(psi1*psi2-psi2-psi1+1)*chi^3*zeta^6-2*chi^3*zeta^4-3*chi^3*zeta^2+...
    (psi1+psi2-3*psi1*psi2+1)*chi^2*zeta^4+2*chi^2*zeta^2+chi^2+3*psi1*psi2*chi*zeta^2-psi1*psi2;
    E1 = psi2*chi^3*zeta^4-psi2*chi^2*zeta^2+psi1*psi2*chi*zeta^2-psi1*psi2;
    E2 = chi^5*zeta^6-3*chi^4*zeta^4+(psi1-1)*chi^3*zeta^6+2*chi^3*zeta^4+3*chi^3*zeta^2+(-psi1-1)*chi^2*zeta^4-...
    2*chi^2*zeta^2-chi^2;
    B = E1 / E0;
    V = E2 / E0;
    risk = F1^2*B + (Fstar^2+tau^2)*V + Fstar^2;
end

