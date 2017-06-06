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

function [sig, IF] = power_sig_model(k, a, fst, fend, tst, tend, t)
% k    : power
% a    : amplitude
% fst  : starting frequency
% fend : finishing frequency
% tst  : starting time
% tend : finishing time
% t    : signal time array

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