%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BlockWiseKronProd computes the block wise Kronecker product
% cf. details in: D. Nion, "A Tensor Framework for Nonunitary Joint Block
% Diagonalization," in IEEE Transactions on Signal Processing,
% vol. 59, no. 10, pp. 4585-4594, Oct. 2011.
%
% The original code can be obtained from: http://dimitri.nion.free.fr/
% Original Author : Dimitri Nion, 2011
% Modified by     : H.B, Post-Doc for Prof. Boualem Boashash, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mat = BlockWiseKronProd(A,B,partA,partB)
[C1, ~] = size(A);
[C2, ~] = size(B);

indA = [0 cumsum(partA)];
indB = [0 cumsum(partB)];
indMat = [0 cumsum(partA.*partB)];

 Mat = zeros(C1 * C2, sum(partA .* partB));
for bInd = 1:length(partB)
    Mat(:, indMat(bInd) + 1:indMat(bInd+1)) = kron(A(:, indA(bInd)+ ...
    1:indA(bInd+1)), B(:, indB(bInd) + 1:indB(bInd + 1)));
end

end