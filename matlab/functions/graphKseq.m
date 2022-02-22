function [K_seq, condN, rankN] = graphKseq(L_seq, tol)
% GRAPHKSEQ
% An alternative of graphK for sequences of Laplacian matrix.
% This function takes as input the sequence of the graph laplacian matrices
% L_seq in the form of an (n x n x m) matrices for a sequence of m steps.
% eigenshuffle() is used for evaluating the sequence of eigenvalues and
% eigenvectors. The final result is returned in the form of an (n x n x m)
% matrix K.
[V_seq, D_seq] = eigenshuffle(L_seq);
Vdim = size(V_seq);
W_seq = zeros(Vdim);
condN = zeros([1, Vdim(end)]);
rankN = zeros(size(condN));
K_seq = zeros(Vdim);
for i=1:Vdim(end)
    % Calculating the left eigenvalues first
    Li = L_seq(:,:,i);
    Vi = V_seq(:,:,i); 
    Wi = inv(Vi)';
    Ki = zeros(size(Li));
    condN(i) = cond(Vi);
    W_seq(:,:,i) = Wi;
    ri = rank(Li);
    rankN(i) = ri;
    m = Vdim(1) - ri;
    % Now calculating the outer products of the first m eigenvectors for
    % summation.
    Di = D_seq(:,i);
    msk = (abs(Di) < tol);
    Wi_z = Wi(:,msk);
    Vi_z = Vi(:,msk);
    for j = 1:m
        wj = Wi_z(:,j);
        vj = Vi_z(:,j);
        Ki = Ki + (vj * wj');
    end
    K_seq(:,:,i) = Ki;
end
end