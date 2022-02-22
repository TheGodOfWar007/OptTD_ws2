function [P_seq, R_seq] = lyapPseq(L_seq, K_seq, Q, alpha)
%RFROMKSEQ Summary of this function goes here
%   Detailed explanation goes here
Kdim = size(K_seq);
P_seq = zeros(Kdim);
R_seq = zeros(Kdim);
for i=1:Kdim(end)
    Li = L_seq(:,:,i);
    Ki = K_seq(:,:,i);
    Ri = Li + alpha*Ki;
    Pi = lyap(-Ri, Q);
    P_seq(:,:,i) = Pi;
    R_seq(:,:,i) = Ri;
end
end

