function [E] = nrmLDAttack(A, E_NRM, options)
% NRMUATTACK: generates the DoS link Disabling attack matrix corresponding to A having a
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
assert(E_NRM < a_nrm, "The specified attack norm exceeds the adjacency matrix norm.")

Adim = size(A);
n = Adim(1);
E = zeros(n);

if ~ENABLE_SYMMETRY
    msk = (A > TOL);
    V = A(msk)';
    idx = find(msk)';
    perm = randperm(numel(V));
    V = CStack(num2cell(V(perm)));
    idx = CStack(num2cell(idx(perm)));
    while true
        eij = V.pop();
        ij = idx.pop();
        E(ij) = eij;
        if norm(E, 'fro') > E_NRM
            E(ij) = 0;
            break;
        end
    end
else
    U = triu(A);
    Umsk = (U > TOL);
    V = U(Umsk)';
    idx = find(Umsk)';
    perm = randperm(numel(V));
    V = CStack(num2cell(V(perm)));
    idx = CStack(num2cell(idx(perm)));

    while true
        eij = V.pop();
        ij = idx.pop();
        i = mod(ij, n);
        j = (ij - i)/n + 1;
        ji = (i-1)*n + j;
        E(ij) = eij;
        E(ji) = eij;
        if norm(E, 'fro') > E_NRM
            E(ij) = 0;
            E(ji) = 0;
            break;
        end
    end
end

end

