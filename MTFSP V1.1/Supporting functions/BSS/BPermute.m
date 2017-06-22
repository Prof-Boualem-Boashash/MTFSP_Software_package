%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of permutation matrix for JBD
% cf. details in: D. Nion, "A Tensor Framework for Nonunitary Joint Block
% Diagonalization," in IEEE Transactions on Signal Processing,
% vol. 59, no. 10, pp. 4585-4594, Oct. 2011.
%
% The original code can be obtained from: http://dimitri.nion.free.fr/
% Original Author : Dimitri Nion, 2011
% Modified by     : H.B, Post-Doc for Prof. Boualem Boashash, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = BPermute(PermutMat, DiagBlocksSizes)
R = length(DiagBlocksSizes);
N = sum(DiagBlocksSizes);
V = diag(1 ./ sqrt(diag(PermutMat * PermutMat'))) * abs(PermutMat);
IndN = 1:N;
indexP1 = [];

[DiagBlocksSizes, pL]= sort(DiagBlocksSizes,'ascend');  
DiagBlocksSizesSortedAsc = DiagBlocksSizes;

while (~isempty(DiagBlocksSizes))
    Lr = DiagBlocksSizes(1);
    numLr = length(find(DiagBlocksSizes == Lr));
    DiagBlocksSizes(1:numLr) = [];
    N1 = zeros(size(V, 1), 1);
    for i=1:size(V,1)
        vi = sort(V(i,:),'descend');
        N1(i,1) = norm(vi(1:Lr), 'fro') ^ 2;
    end
    
    [~, roWIndVLr] = sort(N1,'descend');
    roWIndVLr = roWIndVLr(1:numLr * Lr);
    W = V(roWIndVLr, roWIndVLr);
    V(:,roWIndVLr)=[];
    V(roWIndVLr,:)=[];

    IndexL1 = 1:numLr * Lr;
    indexPL1 = [];
    for i=1:numLr
        x1 = W(1,:);
        [~, indx1] = sort(x1,'descend');
        indexPL1 = [indexPL1  IndexL1(indx1(1:Lr))];
        IndexL1(indx1(1:Lr)) = [];
        W(:,indx1(1:Lr)) = [];
        W(indx1(1:Lr), :) = [];
    end

    indexP1 = [indexP1  IndN(roWIndVLr(indexPL1.').')];
    IndN(roWIndVLr(indexPL1.').') = [];
end
   
Nc = [0, cumsum(DiagBlocksSizesSortedAsc)];
L = cell(1,R);
for r = 1:R
    L{r} = indexP1(Nc(r) + 1:Nc(r+1));
end

L(:, pL) = L;
indexP1 = cell2mat(L);
P = eye(N, N);
P = P(:, indexP1);
end