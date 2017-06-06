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

function [S, ftrue] = signal_model(n, f_end, f_start, N, Fs)

%% Input Check
fr_l = length(f_end);
fs_l = length(f_start);
if(n <= 0), error('Number of signals must be a positive integer'); end
if(fr_l ~= fs_l || fr_l ~= n || fs_l ~= n)
    error('Frequency rates and initial values must have the same size as n');
end
if(N <= 0), error('Signal length must be a positive integer'); end
if(Fs <= 0), error('Sampling frequency must be a positive integer'); end

%% Signal Generation
S      = zeros(n, N);
IF     = zeros(n,N);
theta  = zeros(n,N);
ftrue  = zeros(n,N);
S(:,1) = 1;
for nn = 1:n
    f_rate      = (f_end(nn)-f_start(nn))/N;
    IF(nn,1)    = f_start(nn);
    theta(nn,1) = f_start(nn);
    for t = 1:N-1
        IF(nn,t+1)    = f_start(nn) + f_rate*t;
        theta(nn,t+1) = IF(nn,t+1) + theta(nn,t);
        S(nn,t+1)     = exp(1j*2*pi*theta(nn,t+1));
    end
    % Source instantaneous frequency indexes
    ftrue(nn,:) = mod(floor(2*N*(f_start(nn)+((0:N-1)/(N-1))*(f_end(nn)-f_start(nn))))-Fs/2, N)+1;
end
end