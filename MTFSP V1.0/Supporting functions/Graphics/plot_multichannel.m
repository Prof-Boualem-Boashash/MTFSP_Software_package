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
%                        Plotting MultiChannel Data
%
% Syntax : [out, tag] = plot_multichannel(x, space, fs, c, lw, legd, loc);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <INPUTs>
% x     : Multichannel Data (ch_n x N).
% space : Vertical spacing between multichannel signals.
% fs    : Sampling frequency (Hz)
% c     : Line colour. see help of 'plot' (default : 'b').
% lw    : Line width. see help of 'plot' (default : 1).
% legd  : A flag to include the legend or not (default : 0).
% loc   : Legend location (default : 'NorthEastOutside').
% <OUTPUTs>
% out   : Plot handle.
% tag   : Channel labels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <EXAMPLE>
% x = rand(30,512);
% [h, tag] = plot_multichannel(x, 1, 1,'b');
% set(gca,'yticklabel',tag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out, tag] = plot_multichannel(x, space, fs, c, lw, legd, loc)

if(nargin < 4)
    c = 'b'; lw = 1; legd = 0;
elseif(nargin < 5)
    lw = 1; legd = 0;
elseif(nargin < 6)
    legd = 0;
elseif(nargin < 7 && legd == 1)
    loc = 'NorthEastOutside';
end

[ch_n, N] = size(x);
t = 0:1/fs:(N-1)/fs;
tag = cell(1,2+ch_n);
tag{1,1} = '';
tag{1,2} = '';
for i = 1:ch_n
    dat = x(i,:);
    out = plot(t, dat+1*space - space.*(i-1),'color',c,'linewidth',lw); hold on
    if(i < 10)
        tag{1,i+2} = ['CH. 0' num2str(i)];
    else
        tag{1,i+2} = ['CH. ' num2str(i)];
    end
end
set(gca,'yticklabel',{},'Ytick',(-ch_n*space):space:(space+(space/2)),'GridLineStyle','-');
grid on; axis([0 t(end) (-ch_n*space)+space (space+(space/2))+2*space]);
if(legd), legend(tag,'Location',loc);end
end