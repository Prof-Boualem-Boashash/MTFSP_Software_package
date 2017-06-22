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
%                   Causality Analysis using Phase-Lock Value
%
% Syntax : S = PLV(D, typ)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% D        : MTFD Matrix, which is a cell array with the size m x m (m is the
%            number of sensors)
% typ      : Estimation method (default : 'standard'):
%            'standard' : Conventional definition of PLV.
%            'extended' : Extended time-frequency definition of PLV using
%             auto and cross terms of the MTFD matrix.
%
% <OUTPUTs>
% S        : Brain connectivity matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% s1 = chirp(0:255,0.4,255,0.1);
% s2 = chirp(0:255,0.1,255,0.4);
% S = [s1 ; s2];
% A = inst_model(2, 10, [10 30], 150, 150/2);
% X = awgn(A*S, 30);
% D_wvd   = mtfd(X,'wvd',255,512);
% for m = 1:10, D_wvd{m,m} = real(D_wvd{m,m}); end
% S1 = PLV(D_wvd, 'standard');
% S2 = PLV(D_wvd, 'extended');
% figure; imagesc(abs(S1)); axis xy;
% figure; imagesc(abs(S2)); axis xy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Sout = PLV(D, typ)

Sout = [];

%% Main Inputs Checkup
if(nargin == 0), fprintf(2,'ERROR (No input TFD)\n'); return;
elseif(nargin == 1)
    if(isempty(D)),  fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    ch_n = length(D);
    typ = 'standard';
elseif(nargin == 2)
    if(isempty(D)),  fprintf(2,'ERROR (Input TFD is empty)\n'); return; end
    ch_n = length(D);
    if(isempty(typ)), typ = 'standard'; end
else
    fprintf(2,'ERROR (Extra inputs, please check the function help)\n'); return;
end

%% Main
if(strcmp(typ,'standard') || strcmp(typ,'s'))
    Sout = zeros(ch_n,ch_n);
    for a = 1:ch_n
        for b = 1:ch_n
            temp = wrapTo2Pi(angle(D{a,a}.*conj(D{b,b})));
            Sout(a,b) = mean2(exp(1j.*temp));
        end
    end
elseif(strcmp(typ,'extended') || strcmp(typ,'e'))
    S = zeros(ch_n,ch_n);
    for a = 1:ch_n
        for b = 1:ch_n
            temp   = wrapTo2Pi(angle(D{a,b})+pi);
            S(a,b) = mean2(exp(1j.*temp));
        end
    end
    v    = diag(S);
    dn   = (v*v').^(0.5);
    Sout = S./dn;
else
    fprintf(2,'ERROR (The selected method is not Defined)\n'); return;
end
end