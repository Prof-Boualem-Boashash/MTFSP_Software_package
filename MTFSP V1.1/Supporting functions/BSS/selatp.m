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
%                Selection of auto and overlap terms
%                   for underdetermined TF BSS
%
% Syntax : [TFauto,TFovlap] = selatp(Dxx_aux,Spd,Olap,Eps,W,Dxx)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% Dxx_aux   : MTFD signals using smoothed TFD.
% Spd       : 'slow' apply denoising; 'fast' without denoising
% Olap      : selecting overlaping points 'auto' or 'supp'
% Eps       : Threshold.
% W         : Whitening matrix
% Dxx       : MTFD signals using MWVD.
%
% <OUTPUTs>
% TFauto    : auto-term MTFD matrices
% TFovlap   : overlap-term MTFD matrices         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TFauto,TFovlap] = selatp(Dxx_aux,Spd,Olap,Eps,W,Dxx)
if nargin < 6
    Dxx = Dxx_aux;
end

% ============================================================
% setup information matrix (TFinfo) for later use
% ============================================================
M = size(Dxx_aux,1); % number of sensors
[T,F] = size(Dxx_aux{1,1}); % numbers of time and frequency samples
TFmatrix = [];
for m1 = 1:M
    temp = [];
    for m2 = 1:M
        temp = [temp Dxx_aux{m1,m2}];
    end
    TFmatrix = [TFmatrix;temp];
end

% The TF Matrix data are stored under a block matrix (Structure format in Matlab)

TFinfo = []; % [t f Eb selection] = [time frequency energy selection]
selection = 1;
overlap = 0;
for t=1:T
    for f = 1:F
        tt = t + (0:M-1)*T;
        ff = f + (0:M-1)*F;
        B = TFmatrix(tt,ff);
        Eb = norm(B,'fro');
        TFinfo = [TFinfo; t f Eb selection overlap];
    end
end

% The TFinfo contains 5 columns, where  each is:
% 1: time indices
% 2: frequency indices
% 3: TF signature for a corrosponding 1 and 2.
% 4: Binary Mask of selected TF points for removing noise.
% 5: Binary Mask of selected TF points for reserving autoterms.

L = size(TFinfo,1); % = TxF

% ============================================================
% remove noise if wanted
% ============================================================
if strcmp(Spd,'slow')
    % do nothing
elseif strcmp(Spd,'fast')
    % choose only TF points in a time-slice with significant energy
    for n = 1:T
        Fn = (n-1)*F+1:n*F;
        cutoff = Eps(1)/100*max(TFinfo(Fn,3)); % maximum energy
        low_Eb_points = find(TFinfo(Fn,3)<cutoff);
        TFinfo((n-1)*T+low_Eb_points,4) = 0; % set selection=0
    end
end

% In this step, we modify the 4th column of TFinfo matrix to indicate if
% the point is kept or not by comparing the energy value to a threshold.
% ============================================================
% remove cross-term points
% ============================================================
TFac = find(TFinfo(:,4)==1); % indices of not-noisy points
LL = length(TFac);

if nargin < 6
    TFmatrix_aux = TFmatrix;
elseif nargin == 6
    TFmatrix_aux = [];
    for m1 = 1:M
        temp = [];
        for m2 = 1:M
            temp = [temp Dxx{m1,m2}];
        end
        TFmatrix_aux = [TFmatrix_aux;temp];
    end
end

for k = 1:LL
    n = TFac(k);
    t = TFinfo(n,1);
    f = TFinfo(n,2);
    tt = t + (0:M-1)*T;
    ff = f + (0:M-1)*F;
    B = TFmatrix_aux(tt,ff);
    if trace(W*B*W')/norm(W*B*W','fro') < Eps(2)
        TFinfo(n,4) = 0;
    end
end

% In this step, we modify the 5th column of TFinfo matrix to indicate if
% the point is kept or not by using auto-term criteria selection
% ============================================================
% now all points represent auto-terms, select overlapped points
% ============================================================
if strcmp(Olap,'auto')
    % do nothing
elseif strcmp(Olap,'supp')
    epsilon = 0.4*M; % seuil empirique
    TFauto = find(TFinfo(:,4)==1);
    for k = 1:length(TFauto)
        n = TFauto(k);
        t = TFinfo(n,1);
        f = TFinfo(n,2);
        tt = t + (0:M-1)*T;
        ff = f + (0:M-1)*F;
        B = TFmatrix_aux(tt,ff);
        tol = norm(B)*epsilon;
        rankB = rank(B,tol); % test du rang par rapport à un seuil
        if rankB > 1
            TFinfo(n,5) = 1;
        end
    end
end
% ============================================================
% finally, write to output the auto-term points and overlapped points
% ============================================================
auto = find(TFinfo(:,4)==1);
TFauto = TFinfo(auto,1:2);
ovlap = find(TFinfo(:,5)==1);
TFovlap = TFinfo(ovlap,1:2);
end
