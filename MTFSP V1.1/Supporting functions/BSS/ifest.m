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
%                         IF Estimation
%
% Syntax : [IF,PP] = ifest(TFD,TFp,order)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% TFD      : MTFD signals
% TFp      : Selected TF points corresponding to auto-terms.
% order    : Degree of the polynomial that fits the data
%
% <OUTPUTs>
% IF       : Instantaneous frequency
% PP       : Fiting polynomial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IF,PP] = ifest(TFD,TFp,order)

[N,~]=size(TFD); % N=signal length
tt=TFp(:,1);
ff=TFp(:,2);

use_peak=0; % 0 or 1, better result
if use_peak % polynomial fit on the peaks
    ttt=[];fff=[];
    for ti=1:N
        ind2=find(tt==ti);
        if ~isempty(ind2)
            energy=[];
            [dummy,fi_max]=max(TFD(ti,ff(ind2)));
            fi=ff(ind2(fi_max));
            ttt=[ttt;ti];fff=[fff;fi];
        end
    end
    tt=ttt;ff=fff;
end
P=polyfit(ff,tt/N*0.5,order);
IF=polyval(P,(0:N-1)');
PP=P';
end
