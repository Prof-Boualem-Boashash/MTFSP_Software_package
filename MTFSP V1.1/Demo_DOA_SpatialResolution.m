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
% Last Modification: 21-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
addpath(genpath('Supporting functions'));

%% Main
ind = 1;
while(ind)
    fprintf('Please select one simulation option from the following:\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('1. Present pre-computed results\n2. Re-compute and present results\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('You can select an option by entering its corrosponding number 1 or 2\n');
    opt = input('');
    if(opt == 1 && isnumeric(opt))
        load('DOA Data\DOA_SpatialResolution.mat');
        %% Plotting
        figure('Color',[1 1 1],'Position',[100, 0, 800 600]);
        plot(ra2-ra1,DOA_music_err,'bo--'); hold on;
        plot(ra2-ra1,DOA_tf_music_err,'bv-.'); hold on;
        plot(ra2-ra1,DOA_esprit_err,'r*-'); hold on;
        plot(ra2-ra1,DOA_tf_esprit_err,'r+-'); grid on;
        legend('TF-MUSIC','MUSIC','ESPRIT','TF-ESPRIT');
        xlabel('\delta \theta'); ylabel('NMSE')
        set(gca,'fontweight','bold','fontsize',12);
        %% Saving
        opt1 = input('Do you want to save the results (Y/N)\n','s');
        if(opt1 == 'y' || opt1 == 'Y')
            print(1,'DOA_MTFD_SaptialResolution','-depsc');
        else
            return
        end
        ind = 0;
    elseif(opt == 2 && isnumeric(opt))
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
        m       = 6;                % Number of Sensors
        ra1     = 10;               % Reception angles of n sources
        ra2     = ra1  +linspace(0.1,7,30);
        SNR     = 25;               % Signal-to-noise Ratio
        lamda   = 150;              % Wavelength
        d_space = lamda/2;          % Element spacing in the array of m antennas
        rng(1);                     % Seed initialization of random noise
        % Detection Parameters
        perc    = 0.4;              % Percentage of the STFD maximum power to select high-energy (t,f) points
        iter_N  = 500;
        theta_N = 3e2;
        theta   = linspace(0, 20, theta_N);
        
        %% Computed variables
        t  = 0:1/fs:(N-1)/fs;       % Time array
        f  = 0:fs/(2*M-1):fs/2;     % Frequency array
        
        %% Signal Generation
        S = signal_model(n, f_end, f_init, N, fs); % n LFM signal generation
        
        %% Channel Mixing Model
        DOA_esprit_err    = zeros(1, length(ra2));
        DOA_music_err     = zeros(1, length(ra2));
        DOA_tf_esprit_err = zeros(1, length(ra2));
        DOA_tf_music_err  = zeros(1, length(ra2));
        for i = 1:length(ra2)
            disp(['Dtheta = ' num2str(ra2(i)), ' deg'])
            %% Mixing
            A = inst_model(n, m, [ra1 ra2(i)], lamda, d_space); % Instantaneous Mixing Model
            %% Initialisation
            DOA_esprit    = zeros(1, n);
            DOA_tf_esprit = zeros(1, n);
            DOA_music     = zeros(1, n);
            DOA_tf_music  = zeros(1, n);
            nmse_music    = zeros(1,n);
            nmse_tf_music = zeros(1,n);
            nmse_esprit   = zeros(1,n);
            nmse_tf_esprit = zeros(1,n);
            for j = 1:iter_N
                X = awgn(A*S, SNR); % Additive White Gaussian Noise for provided SNR
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
                DOA_esprit    = doa(X, n, lamda, d_space,'esprit');
                DOA_tf_esprit = doa(D_s, n, lamda, d_space,'tf-esprit');
                P_music       = doa(X, n, lamda, d_space,'music',theta);
                P_tf_music    = doa(D_s, n, lamda, d_space,'tf-music',theta);
                
                %% MUSIC Peaks Allocation
                [~, loc1] = findpeaks(P_music,'NPeaks',n);
                [~, loc2] = findpeaks(P_tf_music,'NPeaks',n);
                DOA_music = peaks_allocation(loc1, theta, [ra1 ra2(i)]);
                DOA_tf_music = peaks_allocation(loc2, theta, [ra1 ra2(i)]);
                RA = [ra1 ra2(i)];
                for k = 1:n
                    nmse_music(1,k)     = (((DOA_music(k) - RA(k))/RA(k)).^2) + nmse_music(1,k);
                    nmse_tf_music(1,k)  = (((DOA_tf_music(k) - RA(k))/RA(k)).^2) + nmse_tf_music(1,k);
                    nmse_esprit(1,k)    = min([1 (((DOA_esprit(k) - RA(k))/RA(k)).^2)]) + nmse_esprit(1,k);
                    nmse_tf_esprit(1,k) = min([1 (((DOA_tf_esprit(k) - RA(k))/RA(k)).^2)]) + nmse_tf_esprit(1,k);
                end
                disp(100*j/iter_N)
            end
            DOA_music_err(1,i)     = mean(nmse_music);
            DOA_tf_music_err(1,i)  = mean(nmse_tf_music);
            DOA_esprit_err(1,i)    = mean(nmse_esprit);
            DOA_tf_esprit_err(1,i) = mean(nmse_tf_esprit);
        end
        DOA_music_err     = DOA_music_err./iter_N;
        DOA_tf_music_err  = DOA_tf_music_err./iter_N;
        DOA_esprit_err    = DOA_esprit_err./iter_N;
        DOA_tf_esprit_err = DOA_tf_esprit_err./iter_N;
        disp('Finished')
        
        %% saving
        save DOA_SpatialResolution DOA_music_err DOA_tf_music_err DOA_esprit_err DOA_tf_esprit_err ra1 ra2 SNR;
        %% Plotting
        figure('Color',[1 1 1],'Position',[100, 0, 800 600]);
        plot(ra2-ra1,DOA_music_err,'bo--'); hold on;
        plot(ra2-ra1,DOA_tf_music_err,'bv-.'); hold on;
        plot(ra2-ra1,DOA_esprit_err,'r*-'); hold on;
        plot(ra2-ra1,DOA_tf_esprit_err,'r+-'); grid on;
        legend('TF-MUSIC','MUSIC','ESPRIT','TF-ESPRIT');
        xlabel('\delta \theta'); ylabel('NMSE')
        set(gca,'fontweight','bold','fontsize',12);
        %% Saving
        opt1 = input('Do you want to save the results (Y/N)\n','s');
        if(opt1 == 'y' || opt1 == 'Y')
            print(1,'DOA_MTFD_SaptialResolution','-depsc');
        else
            return
        end
        ind = 0;
    else
        ind = 1;
        fprintf(2,'Your selection must be either 1 or 2\n');
    end
end