function ppv = formula2(psi1,psi2,lambda,mu1,mustar,F1,Fstar,tau)
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
    %% PPD variance formula
%     term1 = -1i * nu2 * sqrt(lambda_bar*psi1/psi2) * (F1^2/(1-chi*zeta^2)+Fstar^2+tau^2);
%     term2 = -1i * nu1 * (zeta^2/(1-chi*zeta^2)+1) / sqrt(lambda_bar*psi1*psi2);
%     ppv = term1 * (1+term2);
   ppv = F1^2/(1-chi*zeta^2)+Fstar^2+tau^2;
end