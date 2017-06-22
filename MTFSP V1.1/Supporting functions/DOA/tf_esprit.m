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
%               DOA estimation using TF ESPRIT algorithm
%
% Syntax : DOA = tf_esprit(Ds, n, lamda, d)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% Ds       : Selected MTFD Matrix (m x m).
% n        : Number of sources.
% lamda    : Wavelength of the received signal.
% d        : ULA elements spacing.
%
% <OUTPUTs>
% DOA    : Estimated DOA angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DOA = tf_esprit(Ds, n, lamda, d)
[vv,v]   = svd(Ds);                % Eigenvector Decomposition
[~, ind] = sort(diag(v),'descend'); % Sorting the Eigen values
Vs       = vv(:,ind(1:n));          % Subspace filtering

%% ESPRIT Main
V1 = Vs(1:end-1,:);
V2 = Vs(2:end,:);
Eps = V2\V1;
b = eig(Eps);
DOA = acos(-1j*lamda*log(b)/(2*pi*d));

%% Output
DOA = real(DOA')*180/pi;
for nn = 1:n
    if(DOA(nn) > 90), DOA(nn) = 180 - DOA(nn);
end
end