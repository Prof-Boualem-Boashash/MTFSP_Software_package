%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tensMatProdModeN computes the mode-n tensor-matrix product between
% tensor X and matrix A
% Input: 3rd order tensor X, matrix A , and mode n (1,2, or 3)
% Output: mode-n tensor-matrix product prod
% cf. details in: D. Nion, "A Tensor Framework for Nonunitary Joint Block
% Diagonalization," in IEEE Transactions on Signal Processing,
% vol. 59, no. 10, pp. 4585-4594, Oct. 2011.
%
% The original code can be obtained from: http://dimitri.nion.free.fr/
% Original Author : Dimitri Nion, 2011
% Modified by     : H.B, Post-Doc for Prof. Boualem Boashash, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prod = tensMatProdModeN(X, A, n)

if (~(n==1 || n==2 || n==3))
    error('Mode n should be 1, 2, or 3');
end

[L, M, N] = size(X);
[I, L2] = size(A);

if (L2 ~= size(X, n))
    error('Dimension mismatch between A comumns size and dimension n of tensor X');
end

% Perform necessary reshape/permute for 2D matrix product dim match and
% Perform matrix product and inverse reshape/permute
if (n==1)
    Xprime = reshape(X, L, M*N);
    prod = reshape(A * Xprime, I, M, N);
elseif (n==2)
    Xprime = reshape(permute(X, [2 1 3]), M, L*N);
    prod = permute(reshape(A * Xprime, I, L, N), [2 1 3]);
else
    Xprime = reshape(permute(X, [3 1 2]), N, L*M);
    prod = permute(reshape(A * Xprime, I, L, M), [2 3 1]);
end

end

