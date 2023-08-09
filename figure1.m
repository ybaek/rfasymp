clc; clear; close all;
rng(1);

% Model parameters
mu0 = 1/sqrt(2*pi);
mu1 = .5;
mustar = sqrt((pi-2)/(4*pi));
F1 = 1;
Fstar = 0;
tau = 0;
lambda = 1E-3;
psi2 = 3;
psi1grid = psi2/10:psi2/10:psi2*5;
ratios = psi1grid / psi2;

G = numel(psi1grid);
risk_asymp = zeros(1,G);
ppv_asymp = zeros(1,G);
for g = 1:G
    psi1 = psi1grid(g);
    risk_asymp(g) = formula1(psi1,psi2,lambda,mu1,mustar,F1,Fstar,tau);
    ppv_asymp(g) = formula2(psi1,psi2,lambda,mu1,mustar,F1,Fstar,tau);
end

%% Simulated Experiment
d = 100;
n = psi2 * d;
Ngrid = (n/10):(n/10):(n*5);
G = numel(Ngrid);
I = 20;
m = 150;
risk_sim = zeros(G,I);
ppv_sim = zeros(G,I);
alpha = 0;
e = eye(2,d);

beta1 = normrnd(0,1,[d,1]);
beta1 = F1 * beta1 / sqrt(sum(beta1.^2));
for i = 1:I
    X = normrnd(0,1,[n,d]);
    X = X + alpha*e(2,:);
    X = sqrt(d)*(X ./ sqrt(sum(X.^2,2)));
    xnew = normrnd(0,1,[m,d]);
    xnew = xnew + alpha*e(1,:);
    xnew = sqrt(d)*(xnew ./ sqrt(sum(xnew.^2,2)));
    y = X*beta1 + tau*normrnd(0,1,[n,1]);
    fnew = xnew*beta1;
    for g = 1:G
        N = Ngrid(g);
        [risk,ppv] = simulate(y,X,xnew,fnew,N,lambda);
        %[risk,ppv] = simulate_gc(y,X,xnew,fnew,N,lambda,mu0,mu1,mustar);
        risk_sim(g,i) = risk;
        %risk_sim_gc(g,i) = risk;
        ppv_sim(g,i) = ppv;
        %ppv_sim_gc(g,i) = ppv;
    end
end
risk_mean = mean(risk_sim,2);
risk_std  = std(risk_sim,0,2);
ppv_mean = mean(ppv_sim,2);
ppv_std  = std(ppv_sim,0,2);

figure(2)
errorbar(ratios,ppv_mean,ppv_std,'o')
line(ratios,ppv_asymp,'Color','black')
line(ratios,risk_asymp,'Color','red','LineStyle','--')
ylim([0 2])
xlabel('N/n')
ylabel('Expected PPV')
title(strcat("\lambda=",string(lambda)))
legend('S_{RF}^2 (Numerical)','S_{RF}^2 (Theoretical)','R_{RF} (Theoretical)')
fontsize(gcf,scale=1.4)
