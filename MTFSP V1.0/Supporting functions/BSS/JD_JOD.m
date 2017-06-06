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

function V = JD_JOD(aSTFD, xSTFD, N)

% Initialization
L   = min(size(aSTFD));
nm1 = max(size(aSTFD));
nm2 = max(size(xSTFD));
V   = eye(L); % L is the source number
threshold = 1/sqrt(N)/100;
more = 1;
while more, more=0;
    for p = 1:L-1
        for q = p+1:L
            % Givens rotations
            g1 = [aSTFD(p,p:L:nm1)-aSTFD(q,q:L:nm1);aSTFD(p,q:L:nm1)+aSTFD(q,p:L:nm1);1i*(aSTFD(q,p:L:nm1)-aSTFD(p,q:L:nm1))];
            g2 = [aSTFD(p,p:L:nm2)-aSTFD(q,q:L:nm2);aSTFD(p,q:L:nm2)+aSTFD(q,p:L:nm2);1i*(aSTFD(q,p:L:nm2)-aSTFD(p,q:L:nm2)) ];
            [vcp,d] = eig(real(g1*g1'-g2*g2'));
            [~,Ki]  = sort(diag(d)); % sort the Eigen values in ascending order
            angles  = vcp(:,Ki(1)); if(sign(angles(1)) ~= 0), angles = sign(angles(1))*angles; end
            c  = sqrt(0.5 + angles(1)/2);
            sr = 0.5*(angles(2) - 1j*angles(3))/c;
            sc = conj(sr);
            yes = abs(sr) > threshold;
            more = more | yes ;
            if yes,
                colp1 = aSTFD(:,p:L:nm1);
                colq1 = aSTFD(:,q:L:nm1);
                aSTFD(:,p:L:nm1) = c*colp1 + sr*colq1;
                aSTFD(:,q:L:nm1) = c*colq1 - sc*colp1;
                rowp1 = aSTFD(p,:);
                rowq1 = aSTFD(q,:);
                aSTFD(p,:) = c*rowp1 + sc*rowq1;
                aSTFD(q,:) = c*rowq1 - sr*rowp1;
                colp2 = xSTFD(:,p:L:nm2);
                colq2 = xSTFD(:,q:L:nm2);
                xSTFD(:,p:L:nm2) = c*colp2 + sr*colq2;
                xSTFD(:,q:L:nm2) = c*colq2 - sc*colp2;
                rowp2 = xSTFD(p,:);
                rowq2 = xSTFD(q,:);
                xSTFD(p,:) = c*rowp2 + sc*rowq2;
                xSTFD(q,:) = c*rowq2 - sr*rowp2;
                temp = V(:,p);
                V(:,p) = c*V(:,p) + sr*V(:,q);
                V(:,q) = c*V(:,q) - sc*temp;
            end
        end
    end
end
end