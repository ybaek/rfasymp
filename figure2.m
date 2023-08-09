clc; clear; close all;

% Model parameters
mu0 = 1/sqrt(2*pi);
mu1 = .5;
mustar = sqrt((pi-2)/(4*pi));
F1 = 1;
Fstar = 0;
tau = sqrt(5);
psi2 = 3;
psi1grid = psi2/10:psi2/10:psi2*5;
ratios = psi1grid / psi2;

lambda0 = .5;

G = numel(psi1grid);
risk_opt = zeros(1,G);
lambda_opt = zeros(1,G);
options = optimset('Display','iter','TolX',1e-8);
for g = 1:G
    psi1 = psi1grid(g);
    fun1 = @(l) log(formula1(psi2,psi1,exp(l),mu1,mustar,F1,Fstar,tau));
    [lambda,risk] = fminbnd(fun1,-8,2,options);
    lambda_opt(g) = lambda;
    risk_opt(g) = risk;
end
lambda_opt = exp(lambda_opt);
risk_opt = exp(risk_opt);

%
rng(1);
G = numel(psi1grid);
ppv_asymp_opt = zeros(1,G);
for g = 1:G
    psi1 = psi1grid(g);
    lambda1 = lambda_opt(g);
    ppv_asymp_opt(g) = formula2(psi2,psi1,lambda1,mu1,mustar,F1,Fstar,tau);
end

figure(2)
plot(psi1grid,risk_opt ./ (ppv_asymp_opt - tau^2))
xlabel('\psi_2')
ylabel('R / (S^2 - \tau^2)')

% Computing asymptotically optimal values (for psi1 -> infinity)
[lambda_star,omega0,omega1] = lambda_asymp_opt(psi2,F1,Fstar,tau,mu1,mustar);
rho = F1^2 / (Fstar^2+tau^2);
omega_opt = quant_asymp(psi2,lambda_star,mu1,mustar);
risk_wide = (F1^2 + Fstar^2 + tau^2) * (psi2*rho + omega_opt^2) / ...
    ((1+rho)*(psi2-2*omega_opt*psi2+omega_opt^2*psi2-omega_opt^2)) + Fstar^2;
ppv_wide = F1^2 / (1-omega_opt) + Fstar^2;
yline(risk_wide / ppv_wide,'LineStyle','--','Color',[0 0.4470 0.7410])
yline(risk_wide / ppv_wide,'LineStyle','--','Color',[0.8500 0.3250 0.0980])

% If psi2 -> infinity, optimal lambda is always zero
omega_opt = quant_asymp(psi2,0,mu1,mustar);
risk_lsamp = F1^2 * ((omega_opt^3 - omega_opt^2)*mustar^2/mu1^2 + psi2*omega_opt - psi2) / ...
    ((omega_opt-1)*(psi2 - 2*omega_opt*psi2 + omega_opt^2*psi2 - omega_opt^2)) + Fstar^2;
ppv_lsamp = F1^2 / (1-omega_opt) + Fstar^2;
yline(risk_lsamp / ppv_lsamp,'--b')
