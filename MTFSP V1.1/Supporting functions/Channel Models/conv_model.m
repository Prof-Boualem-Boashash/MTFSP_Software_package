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
%                   Instantaneous Mixing Channel
%
% Syntax : [X, A, S] = conv_model(s, Az, L, LL)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% s         : sources signals.
% Az        : mixing filter matrix.
% L         : filter length.
% LL        : window length.
%
% <OUTPUTs>
% X         : mx1 mixed signals matrix.
% A         : mxn convolutive mixing matrix.
% S         : nx1 source signals matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, A, S] = conv_model(s, Az, L, LL)

N = length(s);
[m, n] = size(Az);

%% Constructing the Source Matrix
S = zeros(n*(L+LL), N-LL-L+1);
for jj = 1:n
    St = toeplitz(s(jj,L+LL:-1:1),s(jj,L+LL:end));
    S((jj-1)*(LL+L)+1:jj*(LL+L),:) = St;
end

%% Constructing the Channel Matrix
Ac   = cell(m, n);
A = zeros(m*LL, n*(LL+L));
for ii = 1:m
    for jj = 1:n
        temp = double(coeffs(Az(ii,jj)));
        Aij = toeplitz([temp(1), zeros(1,LL-1)],[temp zeros(1,L+1-length(temp)) zeros(1,LL-1)]);
        Ac{ii,jj} = Aij;
        A((ii-1)*LL+1:ii*LL,(jj-1)*(LL+L)+1:jj*(LL+L)) = Aij;
    end
end
X = A*S;
end