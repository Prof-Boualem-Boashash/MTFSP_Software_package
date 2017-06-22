%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 Boualem Boashash
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Prof. Boualem Boashash        (boualem.boashash@gmail.com)
%          Dr. Abdeldjalil Aissa-El-Bey  (abdeldjalil.aissaelbey@telecom-bretagne.eu)
%          RA: Md.F.A
%
% The following references should be cited whenever this script is used:
% [1] B. Boashash, A. Aissa-El-Bey, Multisensor Time-Frequency Signal Processing:
%     A tutorial review with illustrations in selected application areas, Digital
%     Signal Processing, In Press.
% [2] B. Boashash, A. Aissa-El-Bey, M. F. Al-Sa'd, Multisensor time-frequency
%     signal processing software Matlab package: An analysis tool for multichannel
%     non-stationary data , SoftwareX, In Press.
%
% Last Modification: 25-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            DOA estimation using Conventional MUSIC algorithm
%
% Syntax : P = music(x, n, m, lamda, d, theta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% x        : Input signals (m x N).
% n        : Number of sources.
% m        : Number of sensors.
% lamda    : Wavelength of the received signal.
% d        : ULA elements spacing.
% theta    : Solution space for the 'MUSIC'.
%
% <OUTPUTs>
% P        : Estimated Spectrum for 'MUSIC'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = music(x, n, m, lamda, d, theta)
theta   = theta/180*pi;
N       = length(x); 
theta_N = length(theta);
Rxx     = (1/N)*(x*x');   % Covariance Matrix
[vv,~]  = svd(Rxx);       % Find the eigenvalues and eigenvectors of Rxx
NN      = vv(:,n+1:m);    % Estimate/Selection of noise subspace

%% MUSIC Main
P = zeros(1,theta_N);
for ii = 1:theta_N
    a_theta = exp(-1j*2*pi*(d/lamda)*cos(theta(ii))*(0:m-1));
    P_temp  = conj(a_theta)*(NN*NN')*a_theta.';
    P(ii) = abs(1/P_temp);
end
P = P/max(P);
end
