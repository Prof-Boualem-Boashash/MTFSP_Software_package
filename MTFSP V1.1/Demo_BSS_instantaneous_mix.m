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
fs      = 1;                % Sampling frequency
N       = 256;              % Number of samples
M       = 256;              % Number of frequency bins
% Signal parameters
n       = 3;                % Number of Source Signals
f_init  = [0.0, 0.5, 0.3];  % LFM initial frequencies of n sources
f_end   = [0.5, 0.0, 0.3];  % LFM end frequencies of n sources
% Reception parameters
m       = 6;                % Number of Sensors
ra      = [5, 15, 25];      % Reception angles of n sources
SNR     = 30;               % Signal-to-noise Ratio
lambda  = 150;              % Wavelength
d_space = lambda/2;         % Element spacing in the array of m antennas
% Number of Selected points in Auto and Cross STFDs
Na      = 11;               % Number of Auto-terms matrices used in JD
Nx      = 1;                % Number of Cross-terms matrices used in JAD
th      = 1.8;              % Thereshold for Auto-terms and Cross-terms selection

%% Computed variables
t  = 0:1/fs:(N-1)/fs;       % Time array
f  = 0:fs/(2*M-1):fs/2;     % Frequency array

%% Signal Generation
S = signal_model(n, f_end, f_init, N, fs); % n LFM signal generation
TFD_S = cell(1,n);
for nn = 1:n
    TFD_S{1,nn} = Xwvd(S(nn,:), S(nn,:), N-1, M);
end

%% Channel Mixing Model
rng(1);                   % Seed initialization of random noise
A = inst_model(n, m, ra, lambda, d_space); % Instantaneous Mixing Model
X = awgn(A*S, SNR);       % Additive White Gaussian Noise for provided SNR
TFD_X = cell(1,m);
for mm = 1:m
    TFD_X{1,mm} = Xwvd(X(mm,:), X(mm,:), N-1, M);
end

%% Whitening
W = whitening(X, n, m);
Z = W*X;

%% Source Separation
% Multisensor Time-Frequency Distribution (MTFD)
D = mtfd(Z, 'wvd', N-1, M);
% Selection of Auto and Cross-MTFDs
[aSTFD_all, xSTFD_all] = select_TFD_Instantaneous(D, n, th);
aSTFD = aSTFD_all(:,1:Na*n);
xSTFD = xSTFD_all(:,end+1-Nx*n:end);
% Joint Diagonalization
V = JD_JOD(aSTFD, xSTFD, N);
r = V'*Z;

%% Arranging Estimated Sources
c = zeros(n ,n);
for i = 1:n
    for j = 1:n
        temp = abs(corrcoef(r(i,:), S(j,:)));
        c(i,j) = temp(2);
    end
end
[~, b] = max(c);
R = r(b,:);
TFD_R = cell(1,n);
for nn = 1:n
    TFD_R{1,nn} = Xwvd(R(nn,:), R(nn,:), N-1, M);
end

%% Plotting
tfr = 4;
figure('Color',[1 1 1],'Position',[100, 80, 850 550]);
ha = tight_subplot(2,n,[0.05 0.01],[0.1 0.08],[0.07 0.07]);
for i = 1:n
    axes(ha(i));
    flatwf(f,t(1:tfr:end),abs(TFD_S{1,i}(:,1:tfr:end))','w','k');
    title(['Source Signal ' num2str(i)],'Fontsize',13);
    if i == 1,
        set(gca,'Xtick',[],'fontweight','bold');
        ylabel('Time (s)','Fontsize',12);
    else
        set(gca,'Ytick',[],'Xtick',[]);
    end
    axes(ha(i+n));
    flatwf(f,t(1:tfr:end),abs(TFD_R{1,i}(:,1:tfr:end))','w','k');
    title(['Estimated Source ' num2str(i)],'Fontsize',13);
    if i == 1,
        set(gca,'fontweight','bold');
        ylabel('Time (s)','Fontsize',12);
        xlabel('Frequency (Hz)','Fontsize',12);
    else
        set(gca,'Ytick',[],'fontweight','bold');
        xlabel('Frequency(Hz)','Fontsize',12);
    end
end
figure('Color',[1 1 1],'Position',[100, 80, 500 450]);
ha = tight_subplot(1,1,[0.05 0.01],[0.12 0.08],[0.12 0.12]);
flatwf(f,t(1:tfr:end),abs(TFD_X{1,3}(:,1:tfr:end))','w','k');
ylabel('Time (s)','Fontsize',12); xlabel('Frequency (Hz)','Fontsize',12);
set(gca,'fontweight','bold'); title(['Received Signal on Sensor 3 (SNR = ' num2str(SNR) ' dB)'],'Fontsize',14);

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Inst_BSS_Fig1','-depsc');
    print(2,'Inst_BSS_Fig2','-depsc');
else
    return
end
