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
fs      = 1;                 % Sampling frequency
N       = 256;               % Number of samples
M       = 256;               % Number of frequency bins
% Signal parameters
n       = 3;                 % Number of Source Signals
f_init  = [0.3, 0.4];        % LFM initial frequencies of n-1 sources
f_end   = [0.3, 0.1];        % LFM end frequencies of n-1 sources
fc      = 0.25;              % carrier frequency of the sinusoidal source
fm      = 3.8;               % modulating frequency of the sinusoidal source
fd      = 0.012;             % frequency deviation of the sinusoidal source
a       = 1;                 % Maximum amplitude of the sinusoidal source
tst     = 0;                 % Initial timing of the sinusoidal source
tfi     = N;                 % End timing of the sinusoidal source

%% Computed variables
t       = 0:1/fs:(N-1)/fs;   % Time array
f       = 0:fs/(2*M-1):fs/2; % Frequency array

%% Signal Generation
s1 = sin_sig_model(fc, fd, fm, a, tst, tfi, t);
s2 = signal_model(n-1, f_end, f_init, N, fs); % n-1 LFM signal generation
S = [s1 ; s2];
 
%% Multisensor Time-Frequency Distribution (MTFD)
D_Xwvd   = mtfd(S, 'wvd', N-1, M);
D_Xspwvd = mtfd(S, 'spwvd','hann', 21, 'hann', 31, M);
D_Xckd   = mtfd(S, 'ckd', 1, 0.075, 0.075, M);

%% Plotting
tfr = 5; cnt = 0;
figure('Color',[1 1 1],'Position',[100, 0, 650 680]);
ha = tight_subplot(n,n,[0.05 0.01],[0.1 0.12],[0.1 0.1]);
for i = 1:n
    for j = 1:n
        cnt = cnt + 1;
        axes(ha(cnt)); flatwf(f,t(1:tfr:end),abs(D_Xwvd{i,j}(:,1:tfr:end))','w','k');
        title(sprintf('W_{%s}',['Z_' num2str(i) 'Z_' num2str(j)]))
        if cnt == 1 || cnt == 4,
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
        elseif cnt == 7,
            set(gca,'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
            xlabel('Frequency (Hz)','Fontsize',12);
        elseif cnt == 8 || cnt == 9,
            set(gca,'Ytick',[],'fontweight','bold');
            xlabel('Frequency (Hz)','Fontsize',12);
        else
            set(gca,'Ytick',[],'Xtick',[]);
        end
    end
end
suptitle('Multisensor Time-Frequency Distributions using the WVD')

cnt = 0;
figure('Color',[1 1 1],'Position',[100, 0, 650 680]);
ha = tight_subplot(n,n,[0.05 0.01],[0.1 0.12],[0.1 0.1]);
for i = 1:n
    for j = 1:n
        cnt = cnt + 1;
        axes(ha(cnt)); flatwf(f,t(1:tfr:end),abs(D_Xckd{i,j}(:,1:tfr:end))','w','k');
        title(sprintf('\\rho_{%s}',['Z_' num2str(i) 'Z_' num2str(j)]))
        if cnt == 1 || cnt == 4,
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
        elseif cnt == 7,
            set(gca,'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
            xlabel('Frequency (Hz)','Fontsize',12);
        elseif cnt == 8 || cnt == 9,
            set(gca,'Ytick',[],'fontweight','bold');
            xlabel('Frequency (Hz)','Fontsize',12);
        else
            set(gca,'Ytick',[],'Xtick',[]);
        end
    end
end
suptitle('Multisensor Time-Frequency Distributions using the CKD')

cnt = 0;
figure('Color',[1 1 1],'Position',[100, 0, 650 680]);
ha = tight_subplot(n,n,[0.05 0.01],[0.1 0.12],[0.1 0.1]);
for i = 1:n
    for j = 1:n
        cnt = cnt + 1;
        axes(ha(cnt)); flatwf(f,t(1:tfr:end),abs(D_Xspwvd{i,j}(:,1:tfr:end))','w','k');
        title(sprintf('\\rho_{%s}',['Z_' num2str(i) 'Z_' num2str(j)]))
        if cnt == 1 || cnt == 4,
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
        elseif cnt == 7,
            set(gca,'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
            xlabel('Frequency (Hz)','Fontsize',12);
        elseif cnt == 8 || cnt == 9,
            set(gca,'Ytick',[],'fontweight','bold');
            xlabel('Frequency (Hz)','Fontsize',12);
        else
            set(gca,'Ytick',[],'Xtick',[]);
        end
    end
end
suptitle('Multisensor Time-Frequency Distributions using the SPWVD')

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'MTFD_Fig1','-depsc');
    print(2,'MTFD_Fig2','-depsc');
else
    return
end
