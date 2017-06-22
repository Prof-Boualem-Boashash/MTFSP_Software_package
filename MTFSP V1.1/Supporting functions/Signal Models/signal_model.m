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
%                      Generating n LFM signals
%
% Syntax : [S, ftrue] = signal_model(n, f_end, f_start, N, Fs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% n        : Number of source signals
% f_end    : This is a 1xn array that contains LFM end frequencies.
% f_start  : This is a 1xn array that contains LFM start frequencies.
% N        : Total number of samples
% Fs       : Sampling frequency
%
% <OUTPUTs>
% S        : LFM source signals. This an nxN matrix.
% ftrue    : Instantaneous frequency of each source signal. This an nxN matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% [S, IF] = signal_model(2, [0.4 0.2], [0.1 0.4], 256, 1);
% TFD1 = Xwvd(S(1,:),S(1,:), 255, 256);
% TFD2 = Xwvd(S(2,:),S(2,:), 255, 256);
% figure; 
% subplot(1,2,1); imagesc(0:1/511:0.5,0:255,abs(TFD1.')); axis xy
% subplot(1,2,2); imagesc(0:1/511:0.5,0:255,abs(TFD2.')); axis xy
% figure; 
% subplot(1,2,1); plot(IF(1,:)./512,0:1:255); axis([0 0.5 0 256]);
% subplot(1,2,2); plot(IF(2,:)./512,0:1:255); axis([0 0.5 0 256]);
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