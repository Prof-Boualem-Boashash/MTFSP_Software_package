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
% Last Modification: 27-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc; warning ('off','MATLAB:nargchk:deprecated');
addpath(genpath('Supporting functions'));

%% Parameters that can be changed by the user
% General parameters
fs      = 1;              % Sampling frequency
N       = 256;            % Number of samples
M       = 256;            % Number of frequency bins
% Signal parameters
LFM_f_init  = 0.2;        % LFM initial frequency
LFM_f_end   = 0.45;       % LFM end frequency
QFM_f_init  = 0.45;       % QFM initial frequency
QFM_f_end   = 0.2;        % QFM end frequency
QFM_t_init  = 0;          % QFM initial timing
QFM_t_end   = N;          % QFM end timing
SFM_fc      = 0.1;        % carrier frequency of the sinusoidal source
SFM_fm      = 2;          % modulating frequency of the sinusoidal source
SFM_fd      = 0.012;      % frequency deviation of the sinusoidal source
SFM_a       = 1;          % Maximum amplitude of the sinusoidal source
SFM_tst     = 0;          % Initial timing of the sinusoidal source
SFM_tfi     = N;          % End timing of the sinusoidal source

%% Computed variables
t    = 0:1/fs:(N-1)/fs;   % Time array
f    = 0:fs/(2*M-1):fs/2; % Frequency array

%% Signal Generation
s1 = signal_model(1, LFM_f_end, LFM_f_init, N, fs);
s2 = power_sig_model(3, 1, QFM_f_init, QFM_f_end, QFM_t_init, QFM_t_end, t);
s3 = sin_sig_model(SFM_fc, SFM_fd, SFM_fm, SFM_a, SFM_tst, SFM_tfi, t);
s = s1 + s2 + s3;
TFD_C = Xckd(s, s, 1, 0.08, 0.15, M);

%% Plotting
tfr = 2; 
figure('Color',[1 1 1],'Position',[100, 80, 550 450]);
flatwf(f,t(1:tfr:end),abs(TFD_C(:,1:tfr:end))','w','k');
set(gca,'fontweight','bold'); title('Multicomponent Nonstationary Signal','fontsize',13)
ylabel('Time (s)','Fontsize',12); xlabel('Frequency (Hz)','Fontsize',12);

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'MS_Fig1','-depsc');
else
    return
end
