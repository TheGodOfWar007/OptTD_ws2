function [A] = nrmUMatrix(n, m, NRM, options)
% NRMUMATRIX: generation of matrix having a specific norm with uniformally
% distributed elements b/w the specified limits. If limits are not
% specified, the default is set between 1 and 2.
% @param [in]          n: dimension of square matrix.
% @param [in]          m: number of non-zero elements.
% @param [in]        NRM: norm of the matrix
% @param [in]          p: for calculating the p-norm. (the arg p in norm(A,p))
% @param [opt] SELF_CONN: self connections? T/F. If False, make sure m <=
% n*(n-1), otherwise self connections are kept.
% @param [opt]     w_min: min provisional weight
% @param [opt]     w_max: max provisional weight
% @param [out]         A: output matrix having m non-zero elements with p-norm
% equal to NRM.
% NOTE: WORKS ONLY WITH FROBENIOUS NORM FOR NOW.

arguments
    n double
    m double
    NRM double
    options.p = 'fro'
    options.w_min double = 1
    options.w_max double = 2
    options.TOL double = 1e-6
    options.SELF_CONN logical = false
    options.SYMMETRIC = false
end
w_min = options.w_min;
w_max = options.w_max;
p = options.p;
TOL = options.TOL;
SELF_CONN = options.SELF_CONN;
ENABLE_SYMMETRY = options.SYMMETRIC;

V = w_min + (w_max - w_min)*rand(1,m);
A = zeros(n);
% Placing the values in A randomly.
perm = randperm(n^2);
idx = perm(1:numel(V));
A(idx) = V;
% c_norm = norm(A, p);
% A = A.*(NRM/c_norm);
% We still need to remove self connections. How to do that when matrix is
% strongly connected ? Take the diagonal elements and give their values to
% the ones which are still zero.
if ~SELF_CONN
    assert(m <= n*(n-1),"The number of non-zero links should be less than N*(N-1) for removing self-connections.")
    msk = (A < TOL);
    idx = find(msk);
    % remove diagonal indexes from idx
    d_idx = (1:(n+1):(n^2));
    idx_common = intersect(idx, d_idx);
    idx = setxor(idx, idx_common);
    while numel(idx) < n
        idx = vertcat(idx, d_idx(end));
    end
    % Proceed as usual.
    idx = idx(randperm(numel(idx)));
    A(idx(1:n)) = diag(A);
    A(d_idx) = 0;
end
if ENABLE_SYMMETRY
    A = (A + A')/2;
    % c_norm = norm(A, p);
    % A = A.*(NRM/c_norm);
end
c_norm = norm(A, p);
A = A.*(NRM/c_norm);
end

