function [L] = graph_laplacian(A)
% Returns the graph laplacian for a digraph defined by the adjacency matrix
% A. 
L = diag(sum(A,2)) - A;
end

