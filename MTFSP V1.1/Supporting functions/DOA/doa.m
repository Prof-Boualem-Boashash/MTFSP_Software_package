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
%                   Direction of Arrival Estimation (DOA)
%
% Syntax : P = doa(x, n, lamda, d_space, typ, varargin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% x        : Input signals or selected MTFD Matrix (m x m).
% n        : Number of sources to be estimated.
% lamda    : Wavelength of the received signal (default : 150).
% d_space  : ULA elements spacing (default : lamda/2).
% typ      : DOA estimation method ('MUSIC', 'ESPRIT', 'TF-MUSIC' or 'TF-ESPRIT').
% varargin : Contains theta solution space for the 'MUSIC' and 'TF-MUSIC'
%           (default : linspace(0,40,1e4)).
%
% <OUTPUTs>
% P_DOA    : Estimated Spectrum for 'MUSIC' and 'TF-MUSIC' or DOA angle
%            for 'ESPRIT' and 'TF-ESPRIT'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% s1 = chirp(0:255,0.4,255,0.1);
% s2 = chirp(0:255,0.1,255,0.4);
% S = [s1 ; s2];
% A = inst_model(2, 4, [10 30], 150, 150/2);
% X = awgn(A*S, 30);
% MUSIC  = doa(X, 2, 150, 150/2,'music',linspace(0,50,1e4));
% ESPRIT = doa(X, 2, 150, 150/2, 'ESPRIT');
% plot(linspace(0,50,1e4), MUSIC);
% disp(ESPRIT);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P_DOA = doa(x, n, lamda, d_space, typ, varargin)
P_DOA = [];
%% Main Inputs Checkup
if(nargin == 0), fprintf(2,'ERROR (No Input Signal)\n'); return; end
if(nargin == 1), fprintf(2,'ERROR (Enter the Number of Sources)\n'); return; end
if(isempty(x)),  fprintf(2,'ERROR (No Input Signal)\n'); return; end
if(isempty(n)),  fprintf(2,'ERROR (Enter the Number of Sources)\n'); return; end
check_main_input(nargin);

%% main
if(strcmp(typ,'MUSIC') || strcmp(typ,'music'))
    if(~isempty(varargin))
        if(isempty(varargin{1})), theta = linspace(0,40,1e4);
        else theta = varargin{1}; end;
    end
    P_DOA = music(x, n, m, lamda, d_space, theta);
elseif(strcmp(typ,'TF-MUSIC') || strcmp(typ,'tf-music'))
    if(~isempty(varargin))
        if(isempty(varargin{1})), theta = linspace(0,40,1e4);
        else theta = varargin{1}; end;
    end
    P_DOA = tf_music(x, n, m, lamda, d_space, theta);
elseif(strcmp(typ,'ESPRIT') || strcmp(typ,'esprit'))
    P_DOA = esprit(x, n, lamda, d_space);
    P_DOA = sort(P_DOA);
elseif(strcmp(typ,'TF-ESPRIT') || strcmp(typ,'tf-esprit'))
    P_DOA = tf_esprit(x, n, lamda, d_space);
    P_DOA = sort(P_DOA);
else
    fprintf(2,'ERROR (The selected DOA method is not Defined)\n'); return;
end

%% Checking Main inputs
    function check_main_input(n)
        [xr, xc] = size(x);
        if(xr ~= xc)
            if(xc < xr), x = x'; m = xc;
            elseif(xc > xr), m = xr; end
            if(n == 2)
                [xr, xc] = size(x);
                if(xc < xr), x = x'; m = xc; lamda = 150; d_space = lamda/2; typ = 'music';theta = linspace(0,40,1e4);
                elseif(xc > xr), m = xr; lamda = 150; d_space = lamda/2; typ = 'music';theta = linspace(0,40,1e4);end
            elseif(n == 3)
                [xr, xc] = size(x);
                if(xc < xr), x = x'; m = xc; if(isempty(lamda)), lamda = 150; end; d_space = lamda/2; typ = 'music';theta = linspace(0,40,1e4);
                elseif(xc > xr), m = xr; if(isempty(lamda)), lamda = 150; end; d_space = lamda/2; typ = 'music'; theta = linspace(0,40,1e4);end
            elseif(n == 4)
                [xr, xc] = size(x);
                if(xc < xr), x = x'; m = xc; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; typ = 'music';theta = linspace(0,40,1e4);
                elseif(xc > xr), m = xr; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; typ = 'music';theta = linspace(0,40,1e4);end
            elseif(n == 5)
                [xr, xc] = size(x);
                if(xc < xr), x = x'; m = xc; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'music'; end;theta = linspace(0,40,1e4);
                elseif(xc > xr), m = xr; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'music'; end; theta = linspace(0,40,1e4);end
            elseif(n == 6)
                [xr, xc] = size(x);
                if(xc < xr), x = x'; m = xc; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'music'; end;
                elseif(xc > xr), m = xr; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'music'; end;end
            else
                [xr, xc] = size(x);
                if(xc < xr), x = x'; m = xc; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'music'; end;
                elseif(xc > xr), m = xr; if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'music'; end; end;
            end
        else
            m = xc;
            if(n == 2)
                lamda = 150; d_space = lamda/2; typ = 'tf-music'; theta = linspace(0,40,1e4);
            elseif(n == 3)
                if(isempty(lamda)), lamda = 150; end; d_space = lamda/2; typ = 'tf-music'; theta = linspace(0,40,1e4);
            elseif(n == 4)
                if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; typ = 'tf-music'; theta = linspace(0,40,1e4);
            elseif(n == 5)
                if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'tf-music'; end; theta = linspace(0,40,1e4);
            else
                if(isempty(lamda)), lamda = 150; end; if(isempty(d_space)), d_space = lamda/2; end; if(isempty(typ)), typ = 'tf-music'; end;
            end
        end
    end
end