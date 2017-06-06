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
%               MultiSensor Time-Frequency Distribution (MTFD)
%
% Syntax : D = mtfd(x, typ, varargin);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% x        : Input signals. These can be real or analytic.
% typ      : Type of MTFD 'WVD', 'PWVD', or 'SPWVD' (default : 'wvd').
% varargin : This is a dynamic input argument that contains the parameters
%            of the 'WVD', 'PWVD' or 'SPWVD'. Type 'help Xwvd', 'help Xpwvd',
%            or 'help Xspwvd' for further details on this variable input.
% <OUTPUTs>
% D        : The MTFD Matrix, which is a cell array with the size m x m (m is the
%            number of sensors).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% s1 = chirp(0:255,0.4,255,0.1);
% s2 = chirp(0:255,0.1,255,0.4);
% TFD_wvd   = mtfd([s1;s2],'wvd',255,512);
% TFD_pwvd  = mtfd([s1;s2],'pwvd','rectwin',31,512);
% TFD_spwvd = mtfd([s1;s2],'spwvd','hann',31,'gausswin',31,2,512);
% figure;
% subplot(2,2,1); imagesc(0:1/1023:1/2,0:255,abs(TFD_wvd{1,1}')); axis xy
% subplot(2,2,2); imagesc(0:1/1023:1/2,0:255,abs(TFD_wvd{1,2}')); axis xy
% subplot(2,2,3); imagesc(0:1/1023:1/2,0:255,abs(TFD_wvd{2,1}')); axis xy
% subplot(2,2,4); imagesc(0:1/1023:1/2,0:255,abs(TFD_wvd{2,2}')); axis xy
% figure;
% subplot(2,2,1); imagesc(0:1/1023:1/2,0:255,abs(TFD_pwvd{1,1}')); axis xy
% subplot(2,2,2); imagesc(0:1/1023:1/2,0:255,abs(TFD_pwvd{1,2}')); axis xy
% subplot(2,2,3); imagesc(0:1/1023:1/2,0:255,abs(TFD_pwvd{2,1}')); axis xy
% subplot(2,2,4); imagesc(0:1/1023:1/2,0:255,abs(TFD_pwvd{2,2}')); axis xy
% figure;
% subplot(2,2,1); imagesc(0:1/1023:1/2,0:255,abs(TFD_spwvd{1,1}')); axis xy
% subplot(2,2,2); imagesc(0:1/1023:1/2,0:255,abs(TFD_spwvd{1,2}')); axis xy
% subplot(2,2,3); imagesc(0:1/1023:1/2,0:255,abs(TFD_spwvd{2,1}')); axis xy
% subplot(2,2,4); imagesc(0:1/1023:1/2,0:255,abs(TFD_spwvd{2,2}')); axis xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = mtfd(x, typ, varargin)
D = [];

%% Main Inputs Checkup
if(nargin == 0), fprintf(2,'ERROR (No input signal)\n'); return;
elseif(nargin == 1)
    if(isempty(x)),  fprintf(2,'ERROR (Input signal is empty)\n'); return; end
    [xr, xc] = size(x);
    if(xc < xr), x = x'; m = xc;
    elseif(xc > xr), m = xr; end
    typ = 'wvd';
elseif(nargin >= 2)
    if(isempty(x)),  fprintf(2,'ERROR (Input signal is empty)\n'); return; end
    [xr, xc] = size(x);
    if(xc < xr), x = x'; m = xc; 
    elseif(xc > xr), m = xr; end
    if(isempty(typ)), typ = 'wvd'; end
end

%% MTFD computation
D = cell(m, m);
for i = 1:m;
    for j = 1:m;
        if(strcmp(typ,'WVD') || strcmp(typ,'wvd'))
            if(nargin == 1 || nargin == 2)
                D{i,j} = Xwvd(x(i,:),x(j,:));
            elseif(nargin == 3)
                D{i,j} = Xwvd(x(i,:),x(j,:),varargin{1});
            elseif(nargin == 4)
                D{i,j} = Xwvd(x(i,:),x(j,:),varargin{1},varargin{2});
            end
        elseif(strcmp(typ,'PWVD') || strcmp(typ,'pwvd'))
            if(nargin == 2)
                D{i,j} = Xpwvd(x(i,:),x(j,:));
            elseif(nargin == 3)
                D{i,j} = Xpwvd(x(i,:),x(j,:),varargin{1});
            elseif(nargin == 4)
                D{i,j} = Xpwvd(x(i,:),x(j,:),varargin{1},varargin{2});
            elseif(nargin == 5)
                D{i,j} = Xpwvd(x(i,:),x(j,:),varargin{1},varargin{2},varargin{3});
            elseif(nargin == 6)
                D{i,j} = Xpwvd(x(i,:),x(j,:),varargin{1},varargin{2},varargin{3},varargin{4});
            end
        elseif(strcmp(typ,'SPWVD') || strcmp(typ,'spwvd'))        
            if(nargin == 2)
                D{i,j} = Xspwvd(x(i,:),x(j,:));
            elseif(nargin == 3)
                D{i,j} = Xspwvd(x(i,:),x(j,:),varargin{1});
            elseif(nargin == 4)
                D{i,j} = Xspwvd(x(i,:),x(j,:),varargin{1},varargin{2});
            elseif(nargin == 5)
                D{i,j} = Xspwvd(x(i,:),x(j,:),varargin{1},varargin{2},varargin{3});
            elseif(nargin == 6)
                D{i,j} = Xspwvd(x(i,:),x(j,:),varargin{1},varargin{2},varargin{3},varargin{4});
            elseif(nargin == 7)
                D{i,j} = Xspwvd(x(i,:),x(j,:),varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
            elseif(nargin == 8)
                D{i,j} = Xspwvd(x(i,:),x(j,:),varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
            elseif(nargin == 9)
                D{i,j} = Xspwvd(x(i,:),x(j,:),varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7});
            end             
        else
            fprintf(2,'ERROR (The selected TFD is not Defined)\n'); return;
        end
    end
end
end
