function [lambda,omega0,omega1] = lambda_asymp_opt(psi2,F1,Fstar,tau,mu1,mustar)
omega0 = quant_asymp(psi2,0,mu1,mustar);
omega1 = quant_asymp(psi2,0,F1,sqrt(Fstar^2 + tau^2)); % Determines phase transition
lambda = 0;
if (omega0 < omega1)
    zeta = mu1 / mustar;
    lambda = (zeta^2*psi2 - zeta^2*omega1*psi2 + zeta^2*omega1 + omega1 - omega1^2) / ...
        ((omega1^2 - omega1) * psi2);
    lambda = lambda * mustar^2; % Need re-scaling
end
end