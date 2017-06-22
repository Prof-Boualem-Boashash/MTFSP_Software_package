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
%                   EEG source localization using TF-MUSIC
%
% Syntax : [Loc, P] = EEG_tf_music(Ds, n, m, theta_N, A, orien)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% Ds       : Selected MTFD Matrix (m x m).
% n        : Number of sources.
% m        : Number of sensors.
% theta_N  : Contains theta solution space for the 'MUSIC'.
% A        : Lead field matrix .
% orien    : Dipole orientation 
%
% <OUTPUTs>
% Loc    : Source localization estimates
% P      : Estimated Spectrum for 'MUSIC'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Loc, P] = EEG_tf_music(Ds, n, m, theta_N, A, orien)
[vv,~]  = svd(Ds);          % Find the eigenvalues and eigenvectors of Rxx
VV      = vv(:,n+1:end);    % Estimate/Selection of noise subspace

%% MUSIC Main
P = zeros(n,theta_N);
for k = 1:theta_N-1
    for ii = 1:n
        gains_test = A(:,n*k).*(repmat(orien(1,ii),m,1));
        for jj = 2:n
            gains_test = gains_test + A(:,n*k+jj-1).*(repmat(orien(jj,ii),m,1));
        end
        P(ii,k) = 1/(gains_test'*VV*VV'*gains_test); 
    end
end
Loc = zeros(1,n);
for i = 1:n
    [~,idx] = max(P(i,:));
    Loc(i) = idx;
end
end