function [A, D] = JointBlockDiag(X, DiagBlocksSizes, GDTolerance, MaxGDIter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Joint Block Diagonalization (JBD) of a TFD
% according to the algorithm proposed in "D. Nion, A Tensor Framework
% for Nonunitary Joint Block Diagonalization, in IEEE Transactions
% on Signal Processing, vol. 59, no. 10, pp. 4585-4594, 2011.
% cf. the paper for further details
% 
% The algorithm is based on the Nonlinear Conjugate Gradient  technique 
% with optimal step size applied on the cost function:
%   phi = (1/2) norm(X - D x1 A x2 conj(A) , 'fro')^2
% where xn is the mode-n tensor matrix product
%--------------------------------------------------------------------------
% INPUTS:  
% - X               :  IxIxK tensor, which holds the set of K matrices to block diagonalize
% - DiagBlocksSizes :  =[L1, L2, ..., LR] size of the diagonal blocks (the way the partition is done)
% - GDTolerance     : tolerance to stop the iterations
% - MaxGDIter       : maximum number of iterations
%-------------------------------------------------------------------------------------------------------------------------------------------
% OUTPUTS: 
% - A       :  (IxN) estimate of the block-diagonalizer
% - D       :  (NxNxK) tensor which holds the exactly block-diagonal matrices, estimated in the LS sense from A
%-------------------------------------------------------------------------------------------------------------------------------------------
% The original code can be obtained from: http://dimitri.nion.free.fr/
% Original Author : Dimitri Nion, 2011
% Modified by     : H.B, Post-Doc for Prof. Boualem Boashash, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Main
if (size(X,1) ~= size(X,2))
	error('TFDs should be square matrices');
end
if (MaxGDIter <= 0)
    error('Please provide a positive MaxGDIter');
end

I = size(X, 1); % Size of each square TFD
K = size(X, 3); % number of TFDs
R = length(DiagBlocksSizes); % Number of blocks in the BD matrix D
N = sum(DiagBlocksSizes);
N2 = sum(DiagBlocksSizes.^2);       
Nc = cumsum([0,DiagBlocksSizes]);
N2c = cumsum([0,DiagBlocksSizes.^2]);

% Tranform X into a 2D matrix, each column representing a TFD (IXI rows, K columns)
X2D = reshape(X,size(X,1)*size(X,2),size(X,3));

% Initialization based on the closed form solution for A and D, etc.
%A = ClosedFormSol(X,DiagBlocksSizes,'herm');
[U,S,V] = svd(reshape(X, I^2, K),'econ'); % SVD
NumSlices = 2;
Xprime = reshape(U(:, 1:min(N2,K)) * S(1:min(N2,K), 1:min(N2,K)) * ...
    V(1:NumSlices, 1:min(N2,K))', I, I, NumSlices);
Xp1 = Xprime(:, :, 1);
Xp2 = Xprime(:,:,2);
[u_Xp1,s_Xp1,v_Xp1] = svd(Xp1,0);
u_Xp1 = u_Xp1(:, 1:N);
s_Xp1 = s_Xp1(1:N, 1:N);
v_Xp1 = v_Xp1(:, 1:N);
[u_Xp2,s_Xp2,v_Xp2] = svd(Xp2,0);
u_Xp2 = u_Xp2(:, 1:N);
s_Xp2 = s_Xp2(1:N, 1:N);
v_Xp2 = v_Xp2(:, 1:N);
Xp12 = u_Xp1 * s_Xp1 * v_Xp1' * v_Xp2 /(s_Xp2) * u_Xp2';
[V_Xp12, Diag_Xp12] = eig(Xp12);
Diag_Xp12 = diag(Diag_Xp12);    
[~, indices] = sort(abs(Diag_Xp12),'descend');
Diag_Xp12 = Diag_Xp12(indices);
V_Xp12 = V_Xp12(:,indices);
if (~isreal(Diag_Xp12) && isreal(X))
    rInd = 1;
    T = [0.5 - 1i * 0.5; 0.5 1i * 0.5];
    while (rInd <= I)
        if (~isreal(Diag_Xp12(rInd)))
            V_Xp12(:, rInd:rInd+1) = V_Xp12(:, rInd:rInd+1) * ...
                diag(exp(-1i * sum(angle(V_Xp12(1, rInd:rInd+1))) / 2)) * T;
            rInd = rInd+2;
        else
            rInd = rInd+1;
        end
    end
    V_Xp12 = real(V_Xp12);
end
V_Xp12 = V_Xp12(:, 1:N);
B=pinv(V_Xp12);
Permut = tensMatProdModeN(tensMatProdModeN(X, B, 1), conj(B), 2);
PermutMat = (1/K) * sum(abs(Permut),3);
[PermutEstim] = BPermute(PermutMat, DiagBlocksSizes);
A = V_Xp12 * PermutEstim;
KronMatA = BlockWiseKronProd(conj(A), A, DiagBlocksSizes, DiagBlocksSizes);
Dprime = KronMatA \ X2D;
D=zeros(N,N,K);
for (r = 1:R)
    D(Nc(r) + 1:Nc(r+1), Nc(r) + 1:Nc(r+1), :) = ...
        reshape(Dprime(N2c(r) + 1:N2c(r+1), :), ...
        DiagBlocksSizes(r), DiagBlocksSizes(r), K);
end

% Run the Gradient Descent Algorithm          
Eprime = X2D - KronMatA * Dprime;
phi = norm(Eprime(:),'fro');    
RunDGIter = 1;
DGConvergence = 0; % Criterion to stop GD before MaxGDIter iterations 

% Initialize gradient matrices for A and D and search direction
GDMat = -0.5 * KronMatA' * Eprime;
GAMat = zeros(I,N); 
E=reshape(Eprime,I,I,K);
for k=1:K
    GAMat = GAMat - E(:,:,k)' * A * D(:,:,k) - E(:,:,k)* A * D(:,:,k)';
end
GAMat=0.5*GAMat;
dA=-GAMat;
dD = -GDMat;
[optimStepA, optimStepD] = OptimStepJBD(X2D, A, Dprime, ...
    -GAMat, -GDMat, DiagBlocksSizes);

A = A - optimStepA * GAMat;
Dprime = Dprime - optimStepD*GDMat;


% Main loop of GD
while (~DGConvergence)
    RunDGIter = RunDGIter+1;
    prevPhi = phi; % keep phi old value

   prevGAMat = GAMat;
   prevGDMat = GDMat;
    D=zeros(N,N,K);
    for r = 1:R
        D(Nc(r) + 1:Nc(r+1), Nc(r) + 1:Nc(r+1), :) = ...
            reshape(Dprime(N2c(r) + 1:N2c(r+1), :), ...
            DiagBlocksSizes(r),DiagBlocksSizes(r), K);
    end
    KronMatA = BlockWiseKronProd(conj(A),A,DiagBlocksSizes,DiagBlocksSizes);
    Eprime = X2D - KronMatA * Dprime;  
    E = reshape(Eprime, I, I, K);
    GDMat = -0.5 * KronMatA' * reshape(E, I*I, K);
    GAMat=zeros(I,N);
    for k=1:K
        GAMat = GAMat - E(:,:,k)'*A*D(:,:,k) - E(:,:,k)*A*D(:,:,k)';
    end
    GAMat=0.5*GAMat;
   NGA = norm(GAMat,'fro')^2;
   NGD = norm(GDMat,'fro')^2;
   prevNGA = norm(prevGAMat,'fro')^2;
   prevNGD = norm(prevGDMat,'fro')^2;
   betap = (NGA + NGD) / (prevNGA + prevNGD);
   betap = max([0, betap]);
   if (abs(real(GAMat(:)' * prevGAMat(:)) + real(GDMat(:)' * ...
           prevGDMat(:))) / (NGA+NGD) > 0.1)
      betap=0;
   end
    % Update GD direction and A and D for the next step
   dA = -GAMat     + betap * dA;
   dD = -GDMat + betap * dD;
   [optimStepA,optimStepD]=OptimStepJBD(X2D, A, Dprime, dA, dD, DiagBlocksSizes);
    A =  A + optimStepA * dA;
    Dprime = Dprime + optimStepD * dD;

    % Check convergence
   phi = norm(E(:),'fro');
   if ((abs(phi - prevPhi) / prevPhi < GDTolerance) || ...
           (RunDGIter == MaxGDIter) || (phi < GDTolerance))
       DGConvergence = 1; % Convergence so stop
   end
end
end