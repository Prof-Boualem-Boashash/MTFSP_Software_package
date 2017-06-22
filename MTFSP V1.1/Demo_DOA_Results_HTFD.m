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
%          Dr. Abdeldjalil Aissa-El-Bey  (abdeldjalil.aissaelbey@imt-atlantique.fr)
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
% Last Modification: 25-05-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
addpath(genpath('Supporting functions'));

%% Parameters that can be changed by the user
% General parameters
fs      = 1;             % Sampling frequency
N       = 512;           % Number of samples
M       = 512;           % Number of frequency bins
% Signal parameters
n       = 2;             % Number of Source Signals
f_init  = [0.05, 0.15];  % LFM initial frequencies of n sources
f_end   = [0.15, 0.25];  % LFM end frequencies of n sources
% Reception and Channel parameters
m       = 8;             % Number of Sensors
lambda  = 150;           % Wavelength
d_space = lambda/2;      % Element spacing in the array of m antennas
SNR     = [-10 10];      % Signal-to-noise Ratio
ra      = [-5, 5];       % Reception angles of n sources
% Detection Parameters
perc    = 0.4;           % Percentage of the STFD maximum power to select high-energy (t,f) points
Niter   = 10;            % Number of iterations per SNR
theta_N = 1e4;           % Theta array
rng(4);                  % Seed initialization of random noise

%% Signal Generation
S = signal_model(n, f_end, f_init, N, fs); % n LFM signal generation

%% Channel Mixing Model
A = inst_model(n, m, ra, lambda, d_space); % Instantaneous Mixing Model

%% Main
P_music        = zeros(theta_N, Niter);
P_tf_music_ckd = zeros(theta_N, Niter);
P_tf_music_wvd = zeros(theta_N, Niter);
P_music_avg        = zeros(length(SNR),theta_N);
P_tf_music_ckd_avg = zeros(length(SNR),theta_N);
P_tf_music_wvd_avg = zeros(length(SNR),theta_N);
theta   = linspace(-40, 40, theta_N);
for i = 1:length(SNR)
    for jj = 1:Niter
        X = awgn(A*S, SNR(i));     % Additive White Gaussian Noise for provided SNR
        %% TF-MUSIC Using CKD
        D   = mtfd(X, 'ckd',1, 0.3, 0.3, M);
        %%% Averaged Auto-TFD
        D_avg = zeros(M, N);
        for mm = 1:m, D_avg = D{mm,mm} + D_avg; end
        D_avg = D_avg./m;
        %%% Selection of high-energy (t,f) points
        thr = perc*max(max(D_avg));
        Tr = abs(D_avg) >= thr;
        [F_trace, ~] = find(Tr);
        n_p = length(F_trace);
        D_s = zeros(m, m);
        for m1 = 1:m
            for m2 = 1:m
                D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
            end
        end
        %%% DOA Estimation
        P_tf_music_ckd(:,jj) = tf_music(D_s, n, m, lambda, d_space, theta);
        
        %% MUSIC and TF-MUSIC Using WVD
        D   = mtfd(X, 'wvd', N-1, M);
        %%% Averaged Auto-TFD
        D_avg = zeros(M, N);
        for mm = 1:m, D_avg = D{mm,mm} + D_avg; end
        D_avg = D_avg./m;
        %%% Selection of high-energy (t,f) points
        thr = perc*max(max(D_avg));
        Tr = abs(D_avg) >= thr;
        [F_trace, ~] = find(Tr);
        n_p = length(F_trace);
        D_s = zeros(m, m);
        for m1 = 1:m
            for m2 = 1:m
                D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
            end
        end
        %%% DOA Estimation
        P_music(:,jj)    = music(X, n, m, lambda, d_space, theta);
        P_tf_music_wvd(:,jj) = tf_music(D_s, n, m, lambda, d_space, theta);
    end
    P_music_tmp             = mean(P_music,2);
    P_music_avg(i,:)        = P_music_tmp./max(P_music_tmp);
    P_tf_music_ckd_tmp      = mean(P_tf_music_ckd,2);
    P_tf_music_ckd_avg(i,:) = P_tf_music_ckd_tmp./max(P_tf_music_ckd_tmp);
    P_tf_music_wvd_tmp      = mean(P_tf_music_wvd,2);
    P_tf_music_wvd_avg(i,:) = P_tf_music_wvd_tmp./max(P_tf_music_wvd_tmp);
end

%% Plotting
figure('Color',[1 1 1],'Position',[115, 20, 1100 550]);
ha = tight_subplot(1,length(SNR),[0.05 0.075],[0.12 0.08],[0.07 0.07]);

axes(ha(1));
load('DOA Data\Pmusic_TF-10.mat'); % load results of MDD for -10 dB
plot(theta, 10*log10(P_tf_music_wvd_avg(1,:)),'r-.','LineWidth',3); hold on;
plot(theta, 10*log10(P_tf_music_ckd_avg(1,:)),'g:','LineWidth',3); hold on;
plot(theta, 10*log10(Pmusic_TF(1,:)),'b','LineWidth',2); hold on;
plot(repmat(ra(1),10),linspace(-22,2,10),'k--','linewidth',2); hold on;
plot(repmat(ra(2),10),linspace(-22,2,10),'k--','linewidth',2); hold on;
plot(theta, 10*log10(Pmusic_TF(2,:)),'b','LineWidth',2); grid on;
ylabel('P_M_U_S_I_C (\theta)','Fontsize',16);
axis([-15 15 -5 0.2]); set(gca,'fontweight','bold','fontsize',13);
xlabel('\theta (deg)','Fontsize',16); title('SNR = -10 dB','fontsize',20)
legend('TF-MUSIC Averaged Spectrum based on WVD',...
    'TF-MUSIC Averaged Spectrum based on CKD',...
    'TF-MUSIC Averaged Spectrum based on MDD',...
    'True Angles','Location','Southwest')

axes(ha(2));
load('DOA Data\Pmusic_TF10.mat'); % load results of MDD for 10 dB
plot(theta, 10*log10(P_tf_music_wvd_avg(2,:)),'r-.','LineWidth',3); hold on;
plot(theta, 10*log10(P_tf_music_ckd_avg(2,:)),'g:','LineWidth',3); hold on;
plot(theta, 10*log10(Pmusic_TF(1,:)),'b','LineWidth',2); hold on;
plot(repmat(ra(1),10),linspace(-22,2,10),'k--','linewidth',2); hold on;
plot(repmat(ra(2),10),linspace(-22,2,10),'k--','linewidth',2); hold on;
plot(theta, 10*log10(Pmusic_TF(2,:)),'b','LineWidth',2); grid on;
ylabel('P_M_U_S_I_C (\theta)','Fontsize',16);
axis([-15 15 -20 0.5]); set(gca,'fontweight','bold','fontsize',13);
xlabel('\theta (deg)','Fontsize',16); title('SNR = 10 dB','fontsize',20)

set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'DOA_HTFD','-depsc');
else
    return
end