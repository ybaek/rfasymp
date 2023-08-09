function [risk,ppv,intre,projt,Theta] = simulate_gc(y,X,xnew,fnew,N,lambda,mu0,mu1,mustar)
%SIMULATE Summary of this function goes here
%   Detailed explanation goes here
    n = size(X,1);
    d = size(X,2);
    m = size(xnew,1);
    Theta = normrnd(0,1,[d,N]);
    Theta = sqrt(d)*(Theta ./ sqrt(sum(Theta.^2)));
    Z = mu0 + mu1 * X*Theta / sqrt(d) + mustar * normrnd(0,1,[n,N]); 
    Znew = mu0 + mu1 * xnew*Theta / sqrt(d) + mustar * normrnd(0,1,[m,N]);
    u = lambda*N*n/d;
    hat_mat = Z*Z' + u*eye(n,n);
    hat_a = (Z'*y - Z'*(hat_mat\(Z*Z'*y)))/u;
    proj_mat = (Znew*Znew' - Znew*Z'*(hat_mat\(Z*Znew')))/u;
    risk = sum((fnew - Znew*hat_a).^2)/m;
    intre = (sum(y.^2) - y'*(Z*hat_a))/n;
    projt = trace(proj_mat)/m;
    ppv = intre*(1+projt);
end

