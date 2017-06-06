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

clear; close all; clc;
addpath('EEG data');
addpath(genpath('Supporting functions'));

%% Loading Data
load EEG.mat
s = data; clear data;
[ch_n, N] = size(s);
load Artifact.mat
art = data; clear data;
load Template.mat

%% Parameters that can be changed by the user
n       = 2;          % Number of sources
m       = ch_n;       % Number of Sensors
th      = 1.0;        % Thereshold for Auto-terms and Cross-terms selection
Na      = 7;          % Number of Auto-terms matrices used in JD
Nx      = 1;          % Number of Cross-terms matrices used in JAD
C       = 0.7;        % Artifact Detection Threshold

%% MultiChannel EEG Corruption
art = 1.5*art;
x = s + art;

%% Whitening
w = whitening(x, n, m);
z = w*x;

%% Source Separation
% Multisensor Time-Frequency Distribution (MTFD)
D = mtfd(z, 'wvd', N-1, N);
% Selection of Auto and Cross-STFDs
[aSTFD_all, xSTFD_all] = select_TFD_Instantaneous(D, n, th);
aSTFD = aSTFD_all(:,1:Na*n);
xSTFD = xSTFD_all(:,end+1-Nx*n:end);
% Joint Diagonalization
V = JD_JOD(aSTFD, xSTFD, N);
r = real(V'*z);

%% Artifact Detection and Removal
% Identifying Artifact source signal
c = zeros(1,n);
for i = 1:n
    c(1,i) = abs(corr(template',r(i,:)','type','Spearman'));
end
[~, ind] = max(c);
a_detect = r(ind,:);
% Detecting Artifacts in all Channels
cc = zeros(1,ch_n);
for i = 1:ch_n
    cc(1,i) = abs(corr(x(i,:)',a_detect'));
end
ind = find(cc > C);
% Artifact Estimation
art_est = zeros(ch_n, N);
alpha = mean(x)*a_detect'/(norm(a_detect)^2);
art_est(ind,:) = repmat(alpha*a_detect,length(ind),1);
s_hat = x - art_est;

%% Plotting
M = 10;
indices = 1:M*M;
i1 = 1 + (0:M-1)*M;
cnt = 1; i2 = zeros(1,M*M-M);
for k = 1:M*M
if(~sum(find(i1 == k)))
    i2(cnt) = k;
    cnt = cnt + 1;
end
end
[~, i1] = find(rem(indices,M) ~= 0);
[~, i2] = find(rem(indices,M) == 0);
channel_ind = -1*ones(1,ch_n);
channel_ind(1,ind) = 1;
channel_ind = fliplr(channel_ind);
sh = 100;

figure('Color',[1 1 1],'Position',[100, 0, 900 700]);
subplot(M,M,i1);
[h2, tag] = plot_multichannel(s, sh, fs,'b');
h1 = plot_multichannel(x, sh, fs,'k'); hold on
set(gca,'fontweight','bold','fontsize',16,'yticklabel',tag); xlabel('Time (s)');
legend([h2 h1],'Original Clean EEG','Corrupted EEG',  'Orientation','horizontal','Location','Northwest');
subplot(M,M,i2);
plot(channel_ind,1:ch_n, 'k','linewidth',2);
axis([-1.5 1.5 0 ch_n+2.5]); ylabel('Artifact Contamination Mask'); 
set(gca,'yticklabel',{},'xticklabel',{'NA','','A'},'YAxisLocation','right','XGrid','on',...
    'Ytick',1:ch_n,'GridLineStyle','-','fontweight','bold','fontsize',16);
suptitle('EEG Contamination with MultiChannel Artifacts','Fontweight','bold','Fontsize',18)

figure('Color',[1 1 1],'Position',[100, 0, 900 700]);
subplot(M,M,i1);
h1 = plot_multichannel(s, sh, fs,'b',2); hold on
[h2, tag] = plot_multichannel(s_hat, sh, fs,'r',1);
set(gca,'fontweight','bold','fontsize',16,'yticklabel',tag); xlabel('Time (s)');
legend([h1 h2],'Original Clean EEG','Estimated EEG ','Orientation','horizontal','Location','Northwest');
subplot(M,M,i2);
plot(channel_ind,1:ch_n, 'k','linewidth',2);
axis([-1.5 1.5 0 ch_n+2.5]); ylabel('Artifact Detection Mask'); 
set(gca,'yticklabel',{},'xticklabel',{'NA','','A'},'YAxisLocation','right','XGrid','on',...
    'Ytick',1:ch_n,'GridLineStyle','-','fontweight','bold','fontsize',16);
suptitle('Artifact Detection and Removal using TF-BSS','Fontweight','bold','Fontsize',18)

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'BSS_EEG_Fig1','-depsc');
    print(2,'BSS_EEG_Fig2','-depsc');
else
    return
end
