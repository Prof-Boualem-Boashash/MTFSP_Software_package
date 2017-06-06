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

function [aSTFD, xSTFD]  = select_TFD_Convolutive(STFD, n, Lb, thr)

%% Initialisation
aSTFD = []; % will contain the auto STFDs
xSTFD = []; % will contain the cross STFDs
m = n*Lb;

%% Convert Cell to Matrix
[M, N] = size(STFD{1});
D = zeros(M, N, m, m);
for ii = 1:m
    for jj = 1:m
        D(:,:,ii,jj) = STFD{ii,jj};
    end
end

%% main
count = 1;
for tp = 1:N
    for fp = 1:M
        Ds(:,:) = D(fp,tp,:,:);
        [Vd,~] = eig(Ds);
        Ex = Vd'*Ds*Vd;
        Pex = zeros(n,1);
        for ii=1:n
            Pex(ii) = norm(Ex((ii-1)*Lb+1:ii*Lb,(ii-1)*Lb+1:ii*Lb),'fro');
        end
        if (max(Pex)/norm(Ex,'fro') > 1-thr) % Selection criterion    
            aSTFD(:,:,count) = Ds;
            count = count + 1;
        end
    end
end
end
