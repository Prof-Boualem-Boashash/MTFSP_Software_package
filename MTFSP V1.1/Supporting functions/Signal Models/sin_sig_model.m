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
%                      Generating n SFM signals
%
% Syntax : [sig, ftrue] = sin_sig_model(fc, fd, fm, a, tst, tend, t)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% fc       : center frequency of the sinusoidal FM. This an 1xn array,
%            where n is th number of source signals.
% fd       : frequency deviation in Hz. This a 1xn array.
% fm       : modulation frequency. This a 1xn array.
% a        : source signals amplitude. This a 1xn array.
% tst      : starting times. This a 1xn array.
% tend     : finish times. This a 1xn array.
% t        : time array. This a 1xN array, where N is the total number of
%            samples.
%
% <OUTPUTs>
% sig      : SFM source signals. This an nxN matrix.
% ftrue    : Instantaneous frequency of each source signal. This an nxN matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% [S, IF] = sin_sig_model([0.15 0.35], [0.01 0.005], [2 4], [1 1], [0 0], [255 255], 0:255);
% TFD1 = stft(S(1,:), 256, hanning(63)');
% TFD2 = stft(S(2,:), 256, hanning(21)');
% figure; 
% subplot(1,2,1); imagesc(0:1/511:0.5,0:255,abs(TFD1.')); axis xy
% subplot(1,2,2); imagesc(0:1/511:0.5,0:255,abs(TFD2.')); axis xy
% figure; 
% subplot(1,2,1); plot(IF(1,:),0:1:255); axis([0 0.5 0 256]);
% subplot(1,2,2); plot(IF(2,:),0:1:255); axis([0 0.5 0 256]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sig, IF] = sin_sig_model(fc, fd, fm, a, tst, tend, t)

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