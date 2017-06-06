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
addpath(genpath('Supporting functions'));

%% Parameters that can be changed by the user
% General parameters
fs      = 1;                % Sampling frequency
N       = 512;              % Number of samples
M       = 512;              % Number of frequency bins
% Signal parameters
n       = 2;                % Number of Source Signals
f_init  = [0.0, 0.5];       % LFM initial frequencies of n sources
f_end   = [0.5, 0.0];       % LFM end frequencies of n sources
% Reception parameters
m       = 8;                % Number of Sensors
ra      = [10, 30];         % Reception angles of n sources
SNR_i   = -10:5:15;         % Signal-to-noise Ratio
lamda   = 150;              % Wavelength
d_space = lamda/2;          % Element spacing in the array of m antennas
rng(1);                     % Seed initialization of random noise
% Detection Parameters
perc    = 0.4;              % Percentage of the STFD maximum power to select high-energy (t,f) points
iter_N  = 3000;
theta_N = 1e4;
theta   = linspace(0, 40, theta_N);

%% Computed variables
t  = 0:1/fs:(N-1)/fs;       % Time array
f  = 0:fs/(2*M-1):fs/2;     % Frequency array

%% Signal Generation
S = signal_model(n, f_end, f_init, N, fs); % n LFM signal generation

%% Channel Mixing Model
A = inst_model(n, m, ra, lamda, d_space); % Instantaneous Mixing Model
for i = 1:length(SNR_i)
    disp(['SNR = ' num2str(SNR_i(i)), ' dB'])
    %% Initialisation
    DOA_esprit    = zeros(iter_N, n);
    DOA_tf_esprit = zeros(iter_N, n);
    P_music       = zeros(iter_N, theta_N);
    DOA_music     = zeros(iter_N, n);
    P_tf_music    = zeros(iter_N, theta_N);
    DOA_tf_music  = zeros(iter_N, n);
    for j = 1:iter_N
        X = awgn(A*S, SNR_i(i)); % Additive White Gaussian Noise for provided SNR
        %% Multisensor Time-Frequency Distribution (MTFD)
        D = mtfd(X, 'wvd', N-1, M);
        %% Averaged Auto-TFD
        D_avg = zeros(M, N);
        for mm = 1:m, D_avg = D{mm,mm} + D_avg; end
        D_avg = D_avg./m;
        %% Selection of high-energy (t,f) points
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
        %% DOA Estimation
        DOA_esprit(j,:)    = doa(X, n, lamda, d_space,'esprit');
        P_music(j,:)       = doa(X, n, lamda, d_space,'music',theta);
        P_tf_music(j,:)    = doa(D_s, n, lamda, d_space,'tf-music',theta);
        DOA_tf_esprit(j,:) = doa(D_s, n, lamda, d_space,'tf-esprit');
        
        %% MUSIC Peaks Allocation
        [~, loc1] = findpeaks(P_music(j,:),'NPeaks',n);
        [~, loc2] = findpeaks(P_tf_music(j,:),'NPeaks',n);
        DOA_music(j,:) = peaks_allocation(loc1, theta, ra);
        DOA_tf_music(j,:) = peaks_allocation(loc2, theta, ra);
        disp(100*j/iter_N)
    end
    P_music_avg    = mean(P_music);
    P_music_avg    = P_music_avg./max(P_music_avg);
    P_tf_music_avg = mean(P_tf_music);
    P_tf_music_avg = P_tf_music_avg./max(P_tf_music_avg);
    
    %% Saving
    disp('Saving ...'); SNR = SNR_i(i);
    Path = ['DOA Data\SNR = ' num2str(SNR_i(i)), ' dB']; mkdir(Path);
    save([Path '/Averaged_Spectrum'],'P_music_avg','P_tf_music_avg','theta','ra','-v7.3');
    save([Path '/DOA'],'DOA_music','DOA_tf_music','DOA_esprit','DOA_tf_esprit','SNR','ra','-v7.3')
end
disp('Finished')
