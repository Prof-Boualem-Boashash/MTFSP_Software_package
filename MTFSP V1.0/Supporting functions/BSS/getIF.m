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

function [IF,sig_order]=getIF(Vmap,TFD,TFp)
% =================================================
% IF estimation
% inputs: obtained classes and their associated TFpoints, TFD values
% outputs: IF, arranging IFs in columns

K=size(Vmap,2);
use_IF_true=0;
if use_IF_true % just for testing
    n = 0:N-1; n=n.';
    f11=0.1; f12=0.05;
    f21=0.3; f22=0.2; 
    f31=0.45; f32=0.35;
    a1=(f12-f11)/N; IF(:,1)=f11+a1*n;
    a2=(f22-f21)/N; IF(:,2)=f21+a2*n;
    a3=(f32-f31)/N; IF(:,3)=f31+a3*n;
else
    for k=1:K
        ind=find(Vmap(:,k)==k);
        IF(:,k)=ifest(TFD{1,k},TFp(ind,:),1);
    end
end

%to have correct permutation, assume sources IFs have their starting
%frequency increased:
%so as IF1(1) < IF2(1) < ... <IFk(1)
need_perm=1; %need permutation
if need_perm
    IFs=IF(1,:); %starting frequencies
    [~,sig_order]=sort(IFs);
    IF=IF(:,sig_order);
end
end
