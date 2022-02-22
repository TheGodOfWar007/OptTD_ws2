function [norm_seq] = getNormSeq(A_seq, p)
% GETNORMSEQ 
% This function takes a square matrix sequence A_seq of the form (n x n x m) and
% returns an array of the sequence of norms of the matrices contained in
% A_seq in order. The length of the returned array will be m.
Adim = size(A_seq);
norm_seq = zeros([1, Adim(end)]);
for i=1:Adim(end)
    Ai = A_seq(:,:,i);
    norm_seq(i) = norm(Ai, p);
end
end

