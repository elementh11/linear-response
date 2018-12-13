function [ H ] = iterativeGreensFunction(H00, H01, eps)
% ------------------------------------------------------------
% iterative method for the Green's function
% J. Phys. F: Met. Phys 15(1985)851-858 M. P. Lopez Sancho, J. M. Lopez Sancho and J. Rubio
% ------------------------------------------------------------
niter = 100;
tol = 1.0e-16;

e0 = H00; e0s = H00; e0ds = H00;
alpha0 = H01; beta0 = H01';

for m = 1:niter
    alpha1 = alpha0 /(eps - e0) * alpha0;
    beta1 = beta0 /(eps - e0) * beta0;
    e1 = e0 + alpha0 /(eps - e0) * beta0 + beta0 /(eps - e0) * alpha0;
    e1s = e0s + alpha1;
    e1ds = e0ds + beta1;
    
    conv = max(max(abs(alpha1)));
    conv2 = max(max(abs(beta1)));
    if(conv<tol)&&(conv2<tol); break,end
    alpha0 = alpha1;
    beta0 = beta1;
    e0 = e1;
    e0s = e1s;
    e0ds = e1ds;
end

H.b = e1; H.b(isnan(H.b)) = 0;
H.s = e1s; H.s(isnan(H.s)) = 0;
H.ds = e1ds; H.ds(isnan(H.ds)) = 0;

