function [K, condN, rankN] = graphK(L, tol)
% GRAPHK
% Calculates the K matrix required for calculating the R matrix for the
% modified strict lyapunov equation for general digraphs. Uses
% eigenshuffle() for sorted eigenvalues and eigenvectors.
[V, D] = eigenshuffle(L);
W = inv(V)';
Vdim = size(V);
condN = cond(V);
rankN = rank(L);
msk = (abs(D) <= tol);
V_z = V(:,msk);
W_z = W(:,msk);
m = Vdim(1) - rankN;
K = zeros('like',L);
for i=1:m
    vi = V_z(:,i);
    wi = W_z(:,i);
    K = K + (vi * wi');
end
end

