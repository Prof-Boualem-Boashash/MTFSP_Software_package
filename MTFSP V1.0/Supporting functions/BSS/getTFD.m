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

function TFD=getTFD(TFp,STFD,win)
% It extracts specific points from an input TFD
K=size(TFp,2);
m=size(STFD,1);
TFD=cell(1,K);
[T,F]=size(STFD{1,1});
for k=1:K
    TFD{1,k} = zeros(T,F);
    L=size(TFp{1,k},1);
    for p=1:L
        t=TFp{1,k}(p,1); f=TFp{1,k}(p,2);
        if win==0
            TFDp=[];
            for m1=1:m
                for m2=1:m
                    TFDp(m1,m2) = STFD{m1,m2}(t,f);
                end
            end                
            TFD{1,k}(t,f) = norm(TFDp); % max eigenvalue
        else % collecting points around selected TF point as well
            we=-win:win; % west->east direction
            ns=-win:win; % north->south direction
            ff=f+we; tt=t+ns;
            ff_w=find(ff<=0); ff(ff_w) = []; 
            ff_e=find(ff>F); ff(ff_e) = [];
            tt_n=find(tt<=0); tt(tt_n) = []; 
            tt_s=find(tt>T); tt(tt_s) = [];
            for t_step = 1:length(tt)
                for f_step=1:length(ff)
                    ti = tt(t_step); fi = ff(f_step);
                    % get the STFD of (tt,ff) point
                    TFDp = [];
                    for m1=1:m
                        for m2=1:m
                            TFDp(m1,m2) = STFD{m1,m2}(ti,fi);
                        end
                    end
                    TFD{1,k}(ti,fi) = norm(TFDp); % max eigenvalue
                end
            end
        end        
    end
end
end
