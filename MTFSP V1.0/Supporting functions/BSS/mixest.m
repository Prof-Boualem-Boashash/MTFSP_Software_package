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

function A=mixest(TFp,STFD,use_stft)
% TFp=cell(1,K)=number of classes
% STFD=cell(m,m); m=number of sensors
%A [mxK]

if nargin==2
    use_stft=0;
end
K=size(TFp,2);
m=size(STFD,1);
A=zeros(m,K);
for k = 1:K
    L=size(TFp{1,k},1); % number of TF points for class K
    D=[]; % contains spatial directions of all vectors in class K
    for p = 1:L
        % get the STFD of (tt,ff) point
        t=TFp{1,k}(p,1);
        f=TFp{1,k}(p,2);
        TFD=0;
        if use_stft % if STFD computed by MWVD
            for mi=1:m
                TFD(mi,mi)=STFD{mi,mi}(t,f);
            end
            Dp=diag(TFD)*exp(-angle(1j*TFD(1,1)));
        else
            for m1=1:m
                for m2=1:m
                    TFD(m1,m2)=STFD{m1,m2}(t,f); % not used for now
                end
            end
            [Dp,~,~]=svd(TFD);
            Dp=Dp.*exp(-1j*angle(Dp(1,1)));
        end
        D = [D Dp];
    end
    Dm=mean(D,2);
    A(:,k)=Dm./Dm(1,1);
end
end