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
%                 MDL test for number of sources estimation
%
% Syntax : NSIG = mdltest(x)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% x        : Input multisensor signal (m x N), m is the number of sensors
%            and N observation length
%
% <OUTPUTs>
% NSIG     : Estimated number of sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NSIG = mdltest(x)
K = size(x,2);
N = size(x,1);
R = 1/K*(x*x');
[~, eigenvalsOut] = eig(R);
eigenvals = real(eigenvalsOut);
eigenvals = sort(diag(eigenvals),'descend');
eigenvals(eigenvals<0) = 0;
Ld = zeros(N,1);
for ii=1:N
    Ld(ii) = K*(N-ii)*log(mean(eigenvals(ii+1:N))/(prod(eigenvals(ii+1:N))^(1/(N-ii)))) + (0.5*(ii*(2*N-ii)+1)*log(K));
end
[~,NSIG] = min(Ld);     
end