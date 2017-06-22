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
% Last Modification: 21-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Instantaneous Mixing Channel
%
% Syntax : A = inst_model(n, m, rec_angle,lambda,d_spacing)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% n         : Number of sources.
% m         : Number of sensors.
% rec_angle : Angle of reception in degrees.
% lambda    : Wavelength of the received signal.
% d_spacing : ULA elements spacing.
%
% <OUTPUTs>
% A         : mxn instantaneous mixing matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = inst_model(n, m, rec_angle,lambda,d_spacing)

%% Input Check
rc_l = length(rec_angle);
if(n <= 0), error('Number of signals must be a positive integer'); end
if(m <= 0), error('Number of sensors must be a positive integer'); end
if(rc_l ~= n), error('Reception angles must have the same size as n'); end

%% Generating A matrix
rec_angle = rec_angle/180*pi;
D = zeros(n, m);    % steering matrix with n row and m column
for k = 1:n
    D(k,:) = exp(-1j*2*pi*(d_spacing/lambda)*cos(rec_angle(k))*(0:m-1));
end
A = D.';

end