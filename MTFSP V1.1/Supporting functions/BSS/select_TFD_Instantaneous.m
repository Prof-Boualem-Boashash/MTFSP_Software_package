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
%                Selection of auto and cross terms
%                    for instantaneous TF BSS
%
% Syntax : [aSTFD, xSTFD]  = select_TFD_Instantaneous(STFD, m, thr)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% STFD     : MTFD signals.
% m        : Number of sensors.
% thr      : Threshold.
%
% <OUTPUTs>
% aSTFD    : auto-term MTFD matrices
% xSTFD    : cross-term MTFD matrices        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [aSTFD, xSTFD]  = select_TFD_Instantaneous(STFD, m, thr)

%% Initialisation
aSTFD = []; % will contain the auto STFDs
xSTFD = []; % will contain the cross STFDs

%% Convert Cell to Matrix
[M, N] = size(STFD{1});
D = zeros(M, N, m, m);
for i = 1:m
    for j = 1:m
        D(:,:,i,j) = abs(STFD{i,j});
    end
end

%% main
for tp = 1:N
    for fp = 1:M
        Ds(:,:) = D(fp,tp,:,:);
        if (trace(Ds)/norm(Ds) > thr) % Selection criterion
            aSTFD = [aSTFD Ds]; % Selection of auto STFDs
        else
            xSTFD = [xSTFD Ds]; % Selection of cross STFDs
        end
    end
end
end
