function [E] = nrmUAttack(A, E_NRM, options)
% NRMUATTACK: generates the DoS link weakening attack matrix corresponding to A having a
% norm equal to E_NRM. The links to be attacked are selected randomly, with
% the condition that the norm of E_NRM stays equal to E_NRM.
% @param [in]     A: the original adjacency matrix.
% @param [in] E_NRM: norm of the attack matrix
% @param [opt in]     p: for calculating the p-norm. (the arg p in norm(A,p))
% @param [opt in]   TOL: tolerance for float comparison. [def: 1e-6]
% @param [out]    E: Link disabling attack matrix having p-norm equal to E_NRM.
% NOTE: WORKS ONLY WITH FROBENIOUS NORM FOR NOW.

arguments
    A 
    E_NRM double
    options.p = 'fro'
    options.TOL double = 1e-6
    options.SYMMETRIC = false
end

p = options.p;
TOL = options.TOL;
ENABLE_SYMMETRY = options.SYMMETRIC;

a_nrm = norm(A, p);
if E_NRM > a_nrm
    E_NRM = a_nrm;
end

% Collect the nz indices
msk = (A > TOL);
V = A(msk)';
idx = find(msk)';
Ve = V.*rand(1, numel(V));
Adim = size(A);
E = zeros(Adim(1));
E(idx) = Ve;
if ENABLE_SYMMETRY
    E = (E + E')/2;
end
c_norm = norm(E, p);
E = E.*(E_NRM/c_norm);
end

