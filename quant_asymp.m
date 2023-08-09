function omega = quant_asymp(psi,lambda,mu1,mustar)
zeta = mu1 / mustar;
lambda_bar = lambda / mustar^2;
a = lambda_bar*psi + 1;
b = psi*zeta^2 - zeta^2 - lambda_bar*psi - 1;
c = psi*zeta^2;
omega = -(b + sqrt(b^2 + 4*a*c)) / (2*a);
end