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
%                       Generating n PFM signals
%
% Syntax : [sig, IF] = power_sig_model(k, a, fst, fend, tst, tend, t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% k        : power law. This an 1xn array,
% a        : source signals amplitude. This a 1xn array.
% tst      : starting times. This a 1xn array.
% tend     : finish times. This a 1xn array.
% t        : time array. This a 1xN array, where N is the total number of
%            samples.
%
% <OUTPUTs>
% sig      : PFM source signals. This an nxN matrix.
% ftrue    : Instantaneous frequency of each source signal. This an nxN matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% [sig, IF] = power_sig_model([3 4], [1 1], [0.05 0.25], [0.25 0.45], [0 0], [255 255], 0:255);
% TFD1 = stft(sig(1,:), 256, hanning(63)');
% TFD2 = stft(sig(2,:), 256, hanning(63)');
% figure; 
% subplot(1,2,1); imagesc(0:1/511:0.5,0:255,abs(TFD1.')); axis xy
% subplot(1,2,2); imagesc(0:1/511:0.5,0:255,abs(TFD2.')); axis xy
% figure; 
% subplot(1,2,1); plot(IF(1,:),0:1:255); axis([0 0.5 0 256]);
% subplot(1,2,2); plot(IF(2,:),0:1:255); axis([0 0.5 0 256]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sig, IF] = power_sig_model(k, a, fst, fend, tst, tend, t)

%% Input Check
N     = length(t);
Nk    = length(k);
Na    = length(a);
Nfst  = length(fst);
Nfend = length(fend);
Ntst  = length(tst);
Ntend = length(tend);
if(~isequal(Nk, Na, Nfst, Nfend, Ntst, Ntend))
    error('Signal parameters must all have the same size')
end

%% main
IF  = zeros(Nk, N);
IP  = zeros(Nk, N);
sig = zeros(Nk, N);
for i = 1:Nk
    IF(i,:) = (fst(i)   + (((fend(i)-fst(i))/(tend(i)^(k(i)-1)))*t.^(k(i)-1))).*(t>=tst(i)&t<=tend(i));
    IP(i,:) = (fst(i)*t + (((fend(i)-fst(i))/(tend(i)^(k(i)-1)))*t.^(k(i)-1)).*t./k(i)).*(t>=tst(i)&t<=tend(i));
    sig(i,:) =(a(i).*exp(1j*2*pi.*IP(i,:)).*(t >= tst(i) & t <= tend(i)));  
end
end