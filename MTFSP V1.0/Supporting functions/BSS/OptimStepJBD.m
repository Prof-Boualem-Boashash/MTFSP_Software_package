%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the optimal step for JBD
% cf. details in: D. Nion, "A Tensor Framework for Nonunitary Joint Block
% Diagonalization," in IEEE Transactions on Signal Processing,
% vol. 59, no. 10, pp. 4585-4594, Oct. 2011.
%
% The original code can be obtained from: http://dimitri.nion.free.fr/
% Original Author : Dimitri Nion, 2011
% Modified by     : H.B, Post-Doc for Prof. Boualem Boashash, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optimStepA, optimStepD]=OptimStepJBD(X, MatA, MatD, ...
    directionA, direcionD, DiagBlocksSizes)

KrondirectionAdirectionA = BlockWiseKronProd(conj(directionA), ...
    directionA, DiagBlocksSizes, DiagBlocksSizes);
KronAdirectionA  = BlockWiseKronProd(conj(MatA), directionA, ...
    DiagBlocksSizes, DiagBlocksSizes);
KrondirectionAA  = BlockWiseKronProd(conj(directionA), MatA, ...
    DiagBlocksSizes, DiagBlocksSizes);
KronAA   = BlockWiseKronProd(conj(MatA), MatA, DiagBlocksSizes, ...
    DiagBlocksSizes);
    
% Auxiliary
alpha1 = KrondirectionAdirectionA*direcionD; alpha1 = alpha1(:);
alpha2 = KrondirectionAdirectionA*MatD; alpha2 = alpha2(:);
alpha = KronAdirectionA+KrondirectionAA;
alpha3 = alpha*direcionD; alpha3 = alpha3(:);
alpha4 = alpha*MatD; alpha4 = alpha4(:);
alpha5 = KronAA*direcionD; alpha5 = alpha5(:);
alpha6 = KronAA*MatD-X; alpha6 = alpha6(:);
beta1 = [alpha2, alpha4, alpha6];
beta2 = [alpha1, alpha3, alpha5];
beta12 = beta2' * beta1;
beta22 = beta2' * beta2;
beta11 = beta1' * beta1;
Poly0 = zeros(1,5);
Poly0(1) = real(beta11(1,1));
Poly0(2) = real(2*beta11(1,2));
Poly0(3) = real(2*beta11(1,3)+beta11(2,2));
Poly0(4) = real(2*beta11(3,2));
Poly0(5) = real(beta11(3,3));
Poly1 = zeros(1,5);
Poly1(1) = real(beta12(1,1));
Poly1(2) = real(beta12(1,2)+beta12(2,1));
Poly1(3) = real(beta12(1,3)+beta12(3,1)+beta12(2,2));
Poly1(4) = real(beta12(3,2)+beta12(2,3));
Poly1(5) = real(beta12(3,3));
Poly2 = zeros(1,5);
Poly2(1) = real(beta22(1,1));
Poly2(2) = real(2*beta22(1,2));
Poly2(3) = real(2*beta22(1,3)+beta22(2,2));
Poly2(4) = real(2*beta22(3,2));
Poly2(5) = real(beta22(3,3));
Poly0Poly2 = conv(Poly0,Poly2);
Poly1Poly1 = conv(real(Poly1),real(Poly1));
polDnum1 = -Poly1Poly1+Poly0Poly2;   % degree 8
polDnum2 = (8:-1:1) .* polDnum1(1:8);
Poly2_1= (4:-1:1) .* Poly2(1:4);
SqRoots = roots(conv(polDnum2, Poly2) - conv(polDnum1, Poly2_1));
SqRoots = SqRoots(find(imag(SqRoots) == 0)); % only reals

Value = polyval(polDnum1,SqRoots)./polyval(Poly2,SqRoots);
optimStepA = SqRoots(find(Value== min(Value), 1));
optimStepA = optimStepA(1);
Value2 = [optimStepA^2 optimStepA 1];
optimStepD = -real((Value2 * beta12 * Value2') / (Value2 * beta22 * Value2'));
end