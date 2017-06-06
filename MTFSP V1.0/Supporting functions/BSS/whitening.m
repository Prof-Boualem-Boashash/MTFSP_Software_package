%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 Boualem Boashash
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
% [2] B. Boashash, A. Aissa-El-Bey, Multisensor time-frequency signal processing
%     software Matlab package: An analysis tool for multichannel non-stationary 
%     data , SoftwareX, In Press.
%
% Last Modification: 02-01-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W = whitening(x, n, m)
%% Input Check
[xr, xc] = size(x);
if(xc < xr), x = x'; end
if(n <= 0), error('Number of signals must be a positive integer'); end
if(m <= 0), error('Number of sensors must be a positive integer'); end

%% Computation of the Whitening Matrix W
R = x*x'/length(x); % Compute the covariance matrix of x
if(m > n) % Overdetermined case and assumes white noise
    % Compute the Eigen decomposition of R
    [U, d] = eig(R);
    [power, k] = sort(diag(real(d))); % sort the Eigen values in ascending order
    sigma = mean(power(1:m-n));
    wl = ones(n,1)./sqrt(power(m-n+1:m)-sigma);
    W = diag(wl)*U(1:m,k(m-n+1:m))'; % Computation of the whitening matrix
else % Underdetermined case and assumes no noise
    W = inv(sqrtm(R)); % Computation of the whitening matrix
end
end