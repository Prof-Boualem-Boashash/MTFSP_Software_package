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

function [sig, IF] = sin_sig_model(fc, fd, fm, a, tst, tend, t)
% fc   : carrier frequency
% fd   : frequency deviation
% fm   : modulation frequency
% a    : amplitude
% tst  : starting time
% tend : finishing time
% t    : signal time array

%% Input Check
N     = length(t);
Nfc   = length(fc);
Nfd   = length(fd);
Nfm   = length(fm);
Ntst  = length(tst);
Ntend = length(tend);
if(~isequal(Nfc, Nfd, Nfm, Ntst, Ntend))
    error('Signal parameters must all have the same size')
end

%% main
fm  = fm/N;
fd  = fd.*N;
IF  = zeros(Nfc, N);
IP  = zeros(Nfc, N);
sig = zeros(Nfc, N);
for i = 1:Nfc
    IF(i,:)  = (fc(i) - (fm(i)/2)*(2*pi)*fd(i).*sin(2*pi*t*(fm(i)/2))).*(t>=tst(i)&t<=tend(i));
    IP(i,:)  = (t.*fc(i) + fd(i).*cos(2*pi*t*(fm(i)/2))).*(t>=tst(i)&t<=tend(i));
    sig(i,:) = a(i).*exp(1j*2*pi*IP(i,:)).*(t>=tst(i)&t<=tend(i));
end
end