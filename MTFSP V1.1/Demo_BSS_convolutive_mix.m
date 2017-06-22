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
fs     = 1;                  % Sampling frequency
N      = 128;                % Number of samples
M      = 256;                % Number of frequency bins
% Signal parameters
n      = 2;                  % Number of Source Signals
f_init = [0.5, 0.0];         % Initial frequencies of n sources
f_end  = [0.0, 0.5];         % End frequencies of n sources
t_init = [0, 0];             % Initial timing of n sources
t_end  = [N, N];             % End timing of n sources
k      = [2, 3];             % Instantaneous Phase order of n sources
a      = [1, 1];             % Maximum amplitude of n sources
% Reception parameters
m      = 4;                  % Number of Sensors
SNR    = 40;                 % Signal-to-noise Ratio
Na     = 8;                  % Number of Auto-terms matrices used in JBD
th     = 0.1;                % Thereshold for Auto-terms and Cross-terms selection

%% Computed variables
t  = 0:1/fs:(N-1)/fs;        % Time array
f  = 0:fs/(2*M-1):fs/2;      % Frequency array

%% Signal Generation
S = power_sig_model(k, a, f_init, f_end, t_init, t_end, t);
TFD_S = cell(1,n);
for nn = 1:n
    TFD_S{1,nn} = Xwvd(S(nn,:), S(nn,:), N-1, M);
end

%% Channel Mixing Model
rng(1);                   % Seed initialization of random noise
L = 1; LL = 2;
syms z1
Az = [1.0 + 0.0*z1  0.85 + 0.10*z1
      0.7 + 0.4*z1  0.25 + 1.00*z1
      1.0 + 0.5*z1  0.70 + 0.85*z1
      0.1 + 0.9*z1  0.35 + 0.4*z1];
[X, A] = conv_model(S, Az, L, LL);
X = awgn(X,SNR);
TFD_X = cell(1,m);
for mm = 1:m*LL
    TFD_X{1,mm} = Xwvd(X(mm,:), X(mm,:), N-LL-1, M);
end

%% Whitening
W = whitening(X, n*(LL+L), m*LL);
Z = W*X;

%% Source Separation
% Multisensor Time-Frequency Distribution (MTFD)
D = mtfd(Z, 'wvd', N-LL-1, M);
% Selection of Auto-STFDs
aSTFD_all = select_TFD_Convolutive(D, n, LL+L, th);
aSTFD =  zeros(n*(LL+L), n*(LL+L),Na);
Mx = floor(size(aSTFD_all,3)/Na);
for ii = 1:Na
    aSTFD(:,:,ii) = mean(aSTFD_all(:,:,(ii-1)*Mx+1:ii*Mx),3);
end
% Joint Block Diagonalization
V = JointBlockDiag(aSTFD, (LL+L)*ones(1,n),1e-9,5000);
r = V'*Z;

%% Arrnaging Estimated signals
TFD_r = cell(1,(n*(LL+L)));
for nn = 1:(n*(LL+L))
    TFD_r{1,nn} = Xwvd(r(nn,:), r(nn,:), N-LL-1, M);
end
c = zeros(n,n*(LL+L));
for i = 1:n
    for j = 1:n*(LL+L)
        c(i,j) = abs(corr(reshape(TFD_S{1,i}(:,1:N-LL),1,M*(N-LL))',reshape(TFD_r{1,j},1,M*(N-LL))'));
    end
end
[~, I] = sort(c,2,'descend');
I = I(:,1:(LL+L));
I = reshape(I,1,6);
R = r(I,:);
TFD_R = cell(1,(n*(LL+L)));
for i = 1:n*(LL+L)
    TFD_R{1,i} = Xwvd(R(i,:), R(i,:), N-LL-1, M);
end

%% Plotting
tfr = 3; cnt = 0;
figure('Color',[1 1 1],'Position',[100, 80, 750 600]);
ha = tight_subplot(n*(LL+L)/2+1,n,[0.05 0.01],[0.1 0.07],[0.08 0.08]);
for i = 1:n*(LL+L+1)
    if i <= n
        axes(ha(i));
        flatwf(f,t(1:tfr:end),abs(TFD_S{1,i}(:,1:tfr:end))','w','k');
        title(['Source Signal ' num2str(i)],'Fontsize',12);
        if i == 1,
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
        else
            set(gca,'Ytick',[],'Xtick',[]);
        end
    else
        axes(ha(i));
        flatwf(f,t(1:tfr:end-LL),abs(TFD_R{1,i-n}(:,1:tfr:end))','w','k');
        if mod(i,2) && i ~= n*(LL+L+1)-1
            cnt = cnt + 1;
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
            title(['Estimated Source 1 v.' num2str(cnt)],'Fontsize',12)
        elseif i == n*(LL+L+1)-1
            cnt = cnt + 1;
            set(gca,'fontweight','bold');
            ylabel('Time (s)','Fontsize',12); xlabel('Frequency (Hz)','Fontsize',12);
            title(['Estimated Source 1 v.' num2str(cnt)],'Fontsize',12)
        elseif i == n*(LL+L+1)
            set(gca,'Ytick',[],'fontweight','bold');
            xlabel('Frequency (Hz)','Fontsize',12);
            title(['Estimated Source 2 v.' num2str(cnt)],'Fontsize',12)
        else
            set(gca,'Ytick',[],'Xtick',[]);
            title(['Estimated Source 2 v.' num2str(cnt)],'Fontsize',12)
        end
    end
end
tfr = 2;
figure('Color',[1 1 1],'Position',[100, 80, 500 450]);
ha = tight_subplot(1,1,[0.05 0.01],[0.12 0.08],[0.12 0.12]);
flatwf(f,t(1:tfr:end-LL),abs(TFD_X{1,2}(:,1:tfr:end))','w','k');
ylabel('Time (s)','Fontsize',12); xlabel('Frequency (Hz)','Fontsize',12);
set(gca,'fontweight','bold'); title(['Received Signal on Sensor 2 (SNR = ' num2str(SNR) ' dB)'],'Fontsize',14);

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Conv_BSS_Fig1','-depsc');
    print(2,'Conv_BSS_Fig2','-depsc');
else
    return
end
