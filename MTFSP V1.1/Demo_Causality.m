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
% Last Modification: 25-04-2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc; warning off
addpath(genpath('Supporting functions'));
addpath(genpath('Lead Field Matrix Data'));
addpath(genpath('Causality Data'));

%% Main
ind = 1;
while(ind)
    fprintf('Please select one option from the following:\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('1. Newborn ');fprintf(2,'(Only Presenting results)\n');
    fprintf('2. Adult\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('You can select an option by entering its corrosponding number 1 or 2\n');
    opt = input('');
    if(opt == 1 && isnumeric(opt))
        ind1 = 1;
        while(ind1)
            fprintf('Please select the signal under analysis:\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('1. Signal 1\n2. Signal 2\n');
            fprintf('-------------------------------------------------------------\n');
            opt1 = input('');
            if(opt1 == 1 && isnumeric(opt1))
                %% Loading
                load('Causality Data\PLV_newborn_sig1.mat');
                %% TFD Generation
                k = 1; j = 1;
                temp1 = Xwvd(sigs(1,:), sigs(1,:), N-1, M);
                temp2 = Xwvd(sigs(2,:), sigs(2,:), N-1, M);
                temp3 = Xwvd(sigs(3,:), sigs(3,:), N-1, M);
                TFD0 = temp1 + temp2 + temp3;
                TFD1 = Xwvd(outputsigs(k,:), outputsigs(j,:), 63, M);
                TFD2 = Xrd(outputsigs(k,:), outputsigs(j,:), 63, 0.15, M);
                TFD3 = Xckd(outputsigs(k,:), outputsigs(j,:), 1, 0.1, 0.1, M);
                %% Plotting
                %%% Thresholding
                perc = 95;  % PLV Thresholding %
                Q    = length(SNR);
                S0_m = cell(1,Q);
                S1_m = cell(1,Q);
                S2_m = cell(1,Q);
                S3_m = cell(1,Q);
                for q = 1:Q
                    S0_m{1,q} = S_stn{1,q}.*(abs(S_stn{1,q}) >= (100-perc)*max(max(abs(S_stn{1,q})))/100);
                    S1_m{1,q} = S_RID{1,q}.*(abs(S_RID{1,q}) >= (100-perc)*max(max(abs(S_RID{1,q})))/100);
                    S2_m{1,q} = S_WVD{1,q}.*(abs(S_WVD{1,q}) >= (100-perc)*max(max(abs(S_WVD{1,q})))/100);
                    S3_m{1,q} = S_CKD{1,q}.*(abs(S_CKD{1,q}) >= (100-perc)*max(max(abs(S_CKD{1,q})))/100);
                end
                %%% TFD Plotting
                t = 0:1/fs:(N-1)/fs;
                f = 0:fs/(2*M-1):fs/2;
                tfr = 2;   % TFR plots resolution
                figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                subplot(2,2,1); flatwf(f,t(1:tfr:end),abs(TFD0(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(a) Ideal Time-Frequency Distribution','fontsize',10);
                subplot(2,2,2); flatwf(f,t(1:tfr:end),abs(TFD2(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(b) Rihaczek Distribution','fontsize',10);
                subplot(2,2,3); flatwf(f,t(1:tfr:end),abs(TFD1(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(c) Wigner-Ville Distribution','fontsize',10);
                subplot(2,2,4); flatwf(f,t(1:tfr:end),abs(TFD3(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(d) Compact Kernel Distribution (CKD)','fontsize',10);
                %%% Brain Maps Plotting
                for q = 1:length(SNR)
                    figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                    subplot(2,2,1), imagesc(abs(S0_m{1,q})); axis xy;
                    title('(a) Standard definition using Rihaczek','fontsize',10); xlabel('Channel'); ylabel('Channel');
                    subplot(2,2,2), imagesc(abs(S1_m{1,q})); axis xy;
                    title('(b) Rihaczek Cross TFD','fontsize',10); xlabel('Channel'); ylabel('Channel');
                    subplot(2,2,3), imagesc(abs(S2_m{1,q})); axis xy;
                    title('(c) Wigner-Ville Distribution','fontsize',10); xlabel('Channel'); ylabel('Channel');
                    subplot(2,2,4), imagesc(abs(S3_m{1,q})); axis xy;
                    title('(d) Compact Kernel Distribution (CKD)','fontsize',10); xlabel('Channel'); ylabel('Channel');
                    if(q ~= 3), caxis([min(min(abs(S3_m{1,q})))+0.05 max(max(abs(S3_m{1,q})))-0.35]); end
                end
                %% Saving Results
                opt5 = input('Do you want to save the results (Y/N)\n','s');
                if(opt5 == 'y' || opt5 == 'Y')
                    print(1,'PLV_TFD_1','-depsc');
                    print(2,'PLV_newborn_sig1_10dB','-depsc');
                    print(3,'PLV_newborn_sig1_30dB','-depsc');
                    print(4,'PLV_newborn_sig1_50dB','-depsc');
                else
                    return
                end
                ind1 = 0;
            elseif(opt1 == 2 && isnumeric(opt1))
                %% Loading
                load('Causality Data\PLV_newborn_sig2.mat');
                %% TFD Generation
                k = 1; j = 1;
                temp1 = Xwvd(sigs(1,:), sigs(1,:), N/4-1, M);
                temp2 = Xwvd(sigs(2,:), sigs(2,:), N/4-1, M);
                temp3 = Xwvd(sigs(3,:), sigs(3,:), N/4-1, M);
                TFD0 = temp1 + temp2 + temp3;
                TFD1 = Xwvd(outputsigs(k,:), outputsigs(j,:), 63, M);
                TFD2 = Xrd(outputsigs(k,:), outputsigs(j,:), 63, 0.15, M);
                TFD3 = Xckd(outputsigs(k,:), outputsigs(j,:), 1, 0.1, 0.1, M);
                %% Plotting
                %%% Thresholding
                perc = 5;  % PLV Thresholding %
                Q = length(SNR);
                S0_m = cell(1,Q);
                S1_m = cell(1,Q);
                S2_m = cell(1,Q);
                S3_m = cell(1,Q);
                for q = 1:Q
                    S0_m{1,q} = S_stn{1,q}.*(abs(S_stn{1,q}) >= (100-perc*2.5)*max(max(abs(S_stn{1,q})))/100);
                    S1_m{1,q} = S_RID{1,q}.*(abs(S_RID{1,q}) >= (100-perc*3)*max(max(abs(S_RID{1,q})))/100);
                    S2_m{1,q} = S_WVD{1,q}.*(abs(S_WVD{1,q}) >= (100-perc)*max(max(abs(S_WVD{1,q})))/100);
                    S3_m{1,q} = S_CKD{1,q}.*(abs(S_CKD{1,q}) >= (100-perc)*max(max(abs(S_CKD{1,q})))/100);
                end
                %%% TFD Plotting
                t = 0:1/fs:(N-1)/fs;
                f = 0:fs/(2*M-1):fs/2;
                tfr = 2;   % TFR plots resolution
                figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                subplot(2,2,1); flatwf(f,t(1:tfr:end),abs(TFD0(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(a) Ideal Time-Frequency Distribution','fontsize',10);
                subplot(2,2,2); flatwf(f,t(1:tfr:end),abs(TFD2(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(b) Rihaczek Distribution','fontsize',10);
                subplot(2,2,3); flatwf(f,t(1:tfr:end),abs(TFD1(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(c) Wigner-Ville Distribution','fontsize',10);
                subplot(2,2,4); flatwf(f,t(1:tfr:end),abs(TFD3(:,1:tfr:end))','w','k');
                xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(d) Compact Kernel Distribution (CKD)','fontsize',10);
                %%% Brain Maps Plotting
                for q = 1:length(SNR)
                    figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                    map = [(1:-0.01:0)',(1:-0.01:0)',(1:-0.01:0)']; colormap(map);
                    subplot(2,2,1), imagesc(abs(S0_m{1,q})); axis xy;
                    title('(a) Standard definition using Rihaczek','fontsize',10); xlabel('Channel'); ylabel('Channel');
                    subplot(2,2,2), imagesc(abs(S1_m{1,q})); axis xy;
                    title('(b) Rihaczek Cross TFD','fontsize',10); xlabel('Channel'); ylabel('Channel');
                    subplot(2,2,3), imagesc(abs(S2_m{1,q})); axis xy;
                    title('(c) Wigner-Ville Distribution','fontsize',10); xlabel('Channel'); ylabel('Channel');
                    subplot(2,2,4), imagesc(abs(S3_m{1,q})); axis xy;
                    title('(d) Compact Kernel Distribution (CKD)','fontsize',10); xlabel('Channel'); ylabel('Channel');
                end
                %% Saving Results
                opt5 = input('Do you want to save the results (Y/N)\n','s');
                if(opt5 == 'y' || opt5 == 'Y')
                    print(1,'PLV_TFD_2','-depsc');
                    print(2,'PLV_newborn_sig2_50dB','-depsc');
                else
                    return
                end
                ind1 = 0;
            else
                ind1 = 1;
                fprintf(2,'Your selection must be either 1 or 2\n');
            end
        end
        ind = 0;
        
    elseif(opt == 2 && isnumeric(opt))
        ind1 = 1;
        while(ind1)
            fprintf('Please select the signal under analysis:\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('1. Signal 1\n2. Signal 2\n');
            fprintf('-------------------------------------------------------------\n');
            opt1 = input('');
            if(opt1 == 1 && isnumeric(opt1))
                ind2 = 1;
                while(ind2)
                    fprintf('Please select one simulation option from the following:\n');
                    fprintf('-------------------------------------------------------------\n');
                    fprintf('1. Present pre-computed results\n2. Re-compute and present results\n');
                    opt2 = input('');
                    if(opt2 == 1 && isnumeric(opt2))
                        %% Loading
                        load('Causality Data\PLV_adult_sig1.mat');
                        %% TFD Generation
                        k = 1; j = 1;
                        temp1 = Xwvd(sigs(1,:), sigs(1,:), N-1, M);
                        temp2 = Xwvd(sigs(2,:), sigs(2,:), N-1, M);
                        temp3 = Xwvd(sigs(3,:), sigs(3,:), N-1, M);
                        TFD0 = temp1 + temp2 + temp3;
                        TFD1 = Xwvd(outputsigs(k,:), outputsigs(j,:), 63, M);
                        TFD2 = Xrd(outputsigs(k,:), outputsigs(j,:), 63, 0.15, M);
                        TFD3 = Xckd(outputsigs(k,:), outputsigs(j,:), 1, 0.1, 0.1, M);
                        %% Plotting
                        %%% Thresholding
                        perc = 95;  % PLV Thresholding %
                        Q    = length(SNR);
                        S0_m = cell(1,Q);
                        S1_m = cell(1,Q);
                        S2_m = cell(1,Q);
                        S3_m = cell(1,Q);
                        for q = 1:Q
                            S0_m{1,q} = S_stn{1,q}.*(abs(S_stn{1,q}) >= (100-perc)*max(max(abs(S_stn{1,q})))/100);
                            S1_m{1,q} = S_RID{1,q}.*(abs(S_RID{1,q}) >= (100-perc)*max(max(abs(S_RID{1,q})))/100);
                            S2_m{1,q} = S_WVD{1,q}.*(abs(S_WVD{1,q}) >= (100-perc)*max(max(abs(S_WVD{1,q})))/100);
                            S3_m{1,q} = S_CKD{1,q}.*(abs(S_CKD{1,q}) >= (100-perc)*max(max(abs(S_CKD{1,q})))/100);
                        end
                        %%% TFD Plotting
                        t = 0:1/fs:(N-1)/fs;
                        f = 0:fs/(2*M-1):fs/2;
                        tfr = 2;   % TFR plots resolution
                        figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                        subplot(2,2,1); flatwf(f,t(1:tfr:end),abs(TFD0(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(a) Ideal Time-Frequency Distribution','fontsize',10);
                        subplot(2,2,2); flatwf(f,t(1:tfr:end),abs(TFD2(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(b) Rihaczek Distribution','fontsize',10);
                        subplot(2,2,3); flatwf(f,t(1:tfr:end),abs(TFD1(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(c) Wigner-Ville Distribution','fontsize',10);
                        subplot(2,2,4); flatwf(f,t(1:tfr:end),abs(TFD3(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(d) Compact Kernel Distribution (CKD)','fontsize',10);
                        %%% Brain Maps Plotting
                        for q = 1:length(SNR)
                            figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                            subplot(2,2,1), imagesc(abs(S0_m{1,q})); axis xy;
                            title('(a) Standard definition using Rihaczek','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,2), imagesc(abs(S1_m{1,q})); axis xy;
                            title('(b) Rihaczek Cross TFD','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,3), imagesc(abs(S2_m{1,q})); axis xy;
                            title('(c) Wigner-Ville Distribution','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,4), imagesc(abs(S3_m{1,q})); axis xy;
                            title('(d) Compact Kernel Distribution (CKD)','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            if(q ~= 3), caxis([min(min(abs(S3_m{1,q})))+0.05 max(max(abs(S3_m{1,q})))-0.35]); end
                        end
                        %% Saving Results
                        opt5 = input('Do you want to save the results (Y/N)\n','s');
                        if(opt5 == 'y' || opt5 == 'Y')
                            print(1,'PLV_TFD_1','-depsc');
                            print(2,'PLV_adult_sig1_10dB','-depsc');
                            print(3,'PLV_adult_sig1_30dB','-depsc');
                            print(4,'PLV_adult_sig1_50dB','-depsc');
                        else
                            return
                        end
                        ind2 = 0;
                    elseif(opt2 == 2 && isnumeric(opt2))
                        %% Parameters that can be changed by the user
                        % Signal Parameters
                        fs  = 20;    % Sampling Frequency
                        T   = 13.8;  % Total Duration
                        % Channel Parameters
                        SNR = [10 30 50];
                        rng(1);
                        % Dipole Source Parameters
                        orientation = [0.5 0.9 0.8; -0.026 -0.36 -0.46 ; -0.24 -0.16 -0.2];
                        locations   = [32836 300 15050];
                        
                        %% Signal Generation
                        t = 1:1/fs:(T*fs-1)/fs;
                        sigs(1,:) = gauspuls(t-T/5,5,0.2);
                        sigs(2,:) = gauspuls(t-T/2,6,0.3);
                        tmp = chirp(t(1:ceil(length(t)*0.2)),1,4,5);
                        sigs(3,ceil(length(t)*0.65):ceil(length(t)*0.65)-1+length(tmp)) = tmp;
                        
                        %% Loading Lead Field Matrix
                        run Combine_Data.m;
                        load LeadfieldMatrix_adult.mat;
                        
                        %% Computed Parameters
                        n      = size(sigs,1);
                        N      = length(sigs);
                        M      = 2^nextpow2(N);
                        ch_n   = size(Gain,1);
                        [A, B] = size(Brain);
                        
                        %% Channel Mixing Model
                        outputsigs = zeros(ch_n, N);
                        for i = 1:n
                            gains = Gain(:,n*locations+i-1).*(repmat(orientation(i,:),ch_n,1));
                            outputsigs = outputsigs + gains*sigs;
                        end
                        clear Gain Brain;
                        
                        %% Main
                        Q = length(SNR);
                        S_stn = cell(1,Q);
                        S_RID = cell(1,Q);
                        S_WVD = cell(1,Q);
                        S_CKD = cell(1,Q);
                        for q = 1:length(SNR)
                            disp(['SNR = ' num2str(SNR(q)) ' dB'])
                            %%% Adding Noise
                            noisyOut = zeros(ch_n, N);
                            for i = 1:ch_n
                                noisyOut(i,:)  = awgn(outputsigs(i,:),SNR(q),'measured');
                            end
                            %%% Cross Channel Causality using standard definition
                            disp('Standard Definition');
                            D = mtfd(noisyOut, 'rd', N/4-1, 0.5, M);
                            for i = 1:ch_n
                                D{i,i} = real(D{i,i});
                            end
                            S_stn{1,q} = PLV(D, 'standard');
                            %%% Cross Channel Causality using Rihaczek Distribution
                            disp('Rihaczek Distribution');
                            S_RID{1,q} = PLV(D, 'extended');
                            clear D;
                            %%% Cross Channel Causality using Wigner-Ville Distribution
                            disp('Wigner-Ville Distribution');
                            D = mtfd(noisyOut, 'wvd', N/4-1, M);
                            for i = 1:ch_n
                                D{i,i} = real(D{i,i});
                            end
                            S_WVD{1,q} = PLV(D, 'extended'); clear D;
                            %% Cross Channel Causality using Compact Support Kernel
                            disp('Compact Support Kernel');
                            D = mtfd(noisyOut, 'ckd', 1, 0.1, 0.1, M);
                            for i = 1:ch_n
                                D{i,i} = real(D{i,i});
                            end
                            S_CKD{1,q} = PLV(D, 'extended'); clear D;
                        end
                        
                        %% Saving Results
                        save('Causality Data\PLV_adult_sig1.mat','S_stn','S_RID','S_WVD','S_CKD','SNR','sigs','outputsigs','N','M','fs');
                        
                        %% Plotting
                        load('Causality Data\PLV_adult_sig1.mat');
                        %% TFD Generation
                        k = 1; j = 1;
                        temp1 = Xwvd(sigs(1,:), sigs(1,:), N-1, M);
                        temp2 = Xwvd(sigs(2,:), sigs(2,:), N-1, M);
                        temp3 = Xwvd(sigs(3,:), sigs(3,:), N-1, M);
                        TFD0 = temp1 + temp2 + temp3;
                        TFD1 = Xwvd(outputsigs(k,:), outputsigs(j,:), 63, M);
                        TFD2 = Xrd(outputsigs(k,:), outputsigs(j,:), 63, 0.15, M);
                        TFD3 = Xckd(outputsigs(k,:), outputsigs(j,:), 1, 0.1, 0.1, M);
                        %% Plotting
                        %%% Thresholding
                        perc = 95;  % PLV Thresholding %
                        Q    = length(SNR);
                        S0_m = cell(1,Q);
                        S1_m = cell(1,Q);
                        S2_m = cell(1,Q);
                        S3_m = cell(1,Q);
                        for q = 1:Q
                            S0_m{1,q} = S_stn{1,q}.*(abs(S_stn{1,q}) >= (100-perc)*max(max(abs(S_stn{1,q})))/100);
                            S1_m{1,q} = S_RID{1,q}.*(abs(S_RID{1,q}) >= (100-perc)*max(max(abs(S_RID{1,q})))/100);
                            S2_m{1,q} = S_WVD{1,q}.*(abs(S_WVD{1,q}) >= (100-perc)*max(max(abs(S_WVD{1,q})))/100);
                            S3_m{1,q} = S_CKD{1,q}.*(abs(S_CKD{1,q}) >= (100-perc)*max(max(abs(S_CKD{1,q})))/100);
                        end
                        %%% TFD Plotting
                        t = 0:1/fs:(N-1)/fs;
                        f = 0:fs/(2*M-1):fs/2;
                        tfr = 2;   % TFR plots resolution
                        figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                        subplot(2,2,1); flatwf(f,t(1:tfr:end),abs(TFD0(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(a) Ideal Time-Frequency Distribution','fontsize',10);
                        subplot(2,2,2); flatwf(f,t(1:tfr:end),abs(TFD2(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(b) Rihaczek Distribution','fontsize',10);
                        subplot(2,2,3); flatwf(f,t(1:tfr:end),abs(TFD1(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(c) Wigner-Ville Distribution','fontsize',10);
                        subplot(2,2,4); flatwf(f,t(1:tfr:end),abs(TFD3(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(d) Compact Kernel Distribution (CKD)','fontsize',10);
                        %%% Brain Maps Plotting
                        for q = 1:length(SNR)
                            figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                            subplot(2,2,1), imagesc(abs(S0_m{1,q})); axis xy;
                            title('(a) Standard definition using Rihaczek','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,2), imagesc(abs(S1_m{1,q})); axis xy;
                            title('(b) Rihaczek Cross TFD','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,3), imagesc(abs(S2_m{1,q})); axis xy;
                            title('(c) Wigner-Ville Distribution','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,4), imagesc(abs(S3_m{1,q})); axis xy;
                            title('(d) Compact Kernel Distribution (CKD)','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            if(q ~= 3), caxis([min(min(abs(S3_m{1,q})))+0.05 max(max(abs(S3_m{1,q})))-0.35]); end
                        end
                        %% Saving Results
                        opt5 = input('Do you want to save the results (Y/N)\n','s');
                        if(opt5 == 'y' || opt5 == 'Y')
                            print(1,'PLV_TFD_1','-depsc');
                            print(2,'PLV_adult_sig1_10dB','-depsc');
                            print(3,'PLV_adult_sig1_30dB','-depsc');
                            print(4,'PLV_adult_sig1_50dB','-depsc');
                        else
                            return
                        end
                        ind2 = 0;
                    else
                        ind2 = 1;
                        fprintf(2,'Your selection must be either 1 or 2\n');
                    end
                end
                ind1 = 0;
            elseif(opt1 == 2 && isnumeric(opt1))
                ind2 = 1;
                while(ind2)
                    fprintf('Please select one simulation option from the following:\n');
                    fprintf('-------------------------------------------------------------\n');
                    fprintf('1. Present pre-computed results\n2. Re-compute and present results\n');
                    opt2 = input('');
                    if(opt2 == 1 && isnumeric(opt2))
                        %% Loading
                        load('Causality Data\PLV_adult_sig2.mat');
                        %% TFD Generation
                        k = 1; j = 1;
                        temp1 = Xwvd(sigs(1,:), sigs(1,:), N/4-1, M);
                        temp2 = Xwvd(sigs(2,:), sigs(2,:), N/4-1, M);
                        temp3 = Xwvd(sigs(3,:), sigs(3,:), N/4-1, M);
                        TFD0 = temp1 + temp2 + temp3;
                        TFD1 = Xwvd(outputsigs(k,:), outputsigs(j,:), 63, M);
                        TFD2 = Xrd(outputsigs(k,:), outputsigs(j,:), 63, 0.15, M);
                        TFD3 = Xckd(outputsigs(k,:), outputsigs(j,:), 1, 0.1, 0.1, M);
                        %% Plotting
                        %%% Thresholding
                        perc = 5;  % PLV Thresholding %
                        Q = length(SNR);
                        S0_m = cell(1,Q);
                        S1_m = cell(1,Q);
                        S2_m = cell(1,Q);
                        S3_m = cell(1,Q);
                        for q = 1:Q
                            S0_m{1,q} = S_stn{1,q}.*(abs(S_stn{1,q}) >= (100-perc*3)*max(max(abs(S_stn{1,q})))/100);
                            S1_m{1,q} = S_RID{1,q}.*(abs(S_RID{1,q}) >= (100-perc*3)*max(max(abs(S_RID{1,q})))/100);
                            S2_m{1,q} = S_WVD{1,q}.*(abs(S_WVD{1,q}) >= (100-perc)*max(max(abs(S_WVD{1,q})))/100);
                            S3_m{1,q} = S_CKD{1,q}.*(abs(S_CKD{1,q}) >= (100-perc)*max(max(abs(S_CKD{1,q})))/100);
                        end
                        %%% TFD Plotting
                        t = 0:1/fs:(N-1)/fs;
                        f = 0:fs/(2*M-1):fs/2;
                        tfr = 2;   % TFR plots resolution
                        figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                        subplot(2,2,1); flatwf(f,t(1:tfr:end),abs(TFD0(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(a) Ideal Time-Frequency Distribution','fontsize',10);
                        subplot(2,2,2); flatwf(f,t(1:tfr:end),abs(TFD2(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(b) Rihaczek Distribution','fontsize',10);
                        subplot(2,2,3); flatwf(f,t(1:tfr:end),abs(TFD1(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(c) Wigner-Ville Distribution','fontsize',10);
                        subplot(2,2,4); flatwf(f,t(1:tfr:end),abs(TFD3(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(d) Compact Kernel Distribution (CKD)','fontsize',10);
                        %%% Brain Maps Plotting
                        for q = 1:length(SNR)
                            figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                            map = [(1:-0.01:0)',(1:-0.01:0)',(1:-0.01:0)']; colormap(map);
                            subplot(2,2,1), imagesc(abs(S0_m{1,q})); axis xy;
                            title('(a) Standard definition using Rihaczek','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,2), imagesc(abs(S1_m{1,q})); axis xy;
                            title('(b) Rihaczek Cross TFD','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,3), imagesc(abs(S2_m{1,q})); axis xy;
                            title('(c) Wigner-Ville Distribution','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,4), imagesc(abs(S3_m{1,q})); axis xy;
                            title('(d) Compact Kernel Distribution (CKD)','fontsize',10); xlabel('Channel'); ylabel('Channel');
                        end
                        %% Saving Results
                        opt5 = input('Do you want to save the results (Y/N)\n','s');
                        if(opt5 == 'y' || opt5 == 'Y')
                            print(1,'PLV_TFD_2','-depsc');
                            print(2,'PLV_adult_sig2_50dB','-depsc');
                        else
                            return
                        end
                        ind2 = 0;
                    elseif(opt2 == 2 && isnumeric(opt2))
                        %% Parameters that can be changed by the user
                        % Channel Parameters
                        SNR = 50;
                        rng(1);
                        % Dipole Source Parameters
                        orientation = [0.5 0.9 0.8; -0.026 -0.36 -0.46 ; -0.24 -0.16 -0.2];
                        locations   = [32836 300 15050];
                        fs          = 1;           % Sampling frequency
                        N           = 256;         % Number of samples
                        M           = 256;         % Number of frequency bins
                        % Signal parameters
                        n           = 3;           % Number of source signals
                        LFM_f_init  = 0.2;         % LFM initial frequency
                        LFM_f_end   = 0.45;        % LFM end frequency
                        QFM_f_init  = [0.45 0.05]; % QFM initial frequency
                        QFM_f_end   = [0.05 0.45]; % QFM end frequency
                        QFM_t_init  = [0 0];       % QFM initial timing
                        QFM_t_end   = [N N];       % QFM end timing
                        
                        %% Computed variables
                        t = 0:1/fs:(N-1)/fs;   % Time array
                        
                        %% Signal Generation
                        s1 = power_sig_model([3 3], [1 1], QFM_f_init, QFM_f_end, QFM_t_init, QFM_t_end, t);
                        s2 = signal_model(1, LFM_f_end, LFM_f_init, N, fs);
                        sigs  = [s1 ; s2];
                        
                        %% Loading Lead Field Matrix
                        run Combine_Data.m;
                        load LeadfieldMatrix_adult.mat;
                        
                        %% Computed Parameters
                        ch_n   = size(Gain,1);
                        [A, B] = size(Brain);
                        
                        %% Channel Mixing Model
                        outputsigs = zeros(ch_n, N);
                        for i = 1:n
                            gains = Gain(:,n*locations+i-1).*(repmat(orientation(i,:),ch_n,1));
                            outputsigs = outputsigs + gains*sigs;
                        end
                        clear Gain Brain;
                        
                        %% Main
                        Q = length(SNR);
                        S_stn = cell(1,Q);
                        S_RID = cell(1,Q);
                        S_WVD = cell(1,Q);
                        S_CKD = cell(1,Q);
                        for q = 1:length(SNR)
                            disp(['SNR = ' num2str(SNR(q)) ' dB'])
                            %%% Adding Noise
                            noisyOut = zeros(ch_n, N);
                            for i = 1:ch_n
                                noisyOut(i,:)  = awgn(outputsigs(i,:),SNR(q),'measured');
                            end
                            %%% Cross Channel Causality using standard definition
                            disp('Standard Definition');
                            D = mtfd(noisyOut, 'rd', N/4-1, 0.5, M);
                            for i = 1:ch_n
                                D{i,i} = real(D{i,i});
                            end
                            S_stn{1,q} = PLV(D, 'standard');
                            %%% Cross Channel Causality using Rihaczek Distribution
                            disp('Rihaczek Distribution');
                            S_RID{1,q} = PLV(D, 'extended');
                            clear D;
                            %%% Cross Channel Causality using Wigner-Ville Distribution
                            disp('Wigner-Ville Distribution');
                            D = mtfd(noisyOut, 'wvd', N/4-1, M);
                            for i = 1:ch_n
                                D{i,i} = real(D{i,i});
                            end
                            S_WVD{1,q} = PLV(D, 'extended'); clear D;
                            %% Cross Channel Causality using Compact Support Kernel
                            disp('Compact Support Kernel');
                            D = mtfd(noisyOut, 'ckd', 1, 0.1, 0.1, M);
                            for i = 1:ch_n
                                D{i,i} = real(D{i,i});
                            end
                            S_CKD{1,q} = PLV(D, 'extended'); clear D;
                        end
                        
                        %% Saving Results
                        save('Causality Data\PLV_adult_sig2.mat','S_stn','S_RID','S_WVD','S_CKD','SNR','sigs','outputsigs','N','M','fs');
                        
                        %% Loading
                        load('Causality Data\PLV_adult_sig2.mat');
                        %% TFD Generation
                        k = 1; j = 1;
                        temp1 = Xwvd(sigs(1,:), sigs(1,:), N/4-1, M);
                        temp2 = Xwvd(sigs(2,:), sigs(2,:), N/4-1, M);
                        temp3 = Xwvd(sigs(3,:), sigs(3,:), N/4-1, M);
                        TFD0 = temp1 + temp2 + temp3;
                        TFD1 = Xwvd(outputsigs(k,:), outputsigs(j,:), 63, M);
                        TFD2 = Xrd(outputsigs(k,:), outputsigs(j,:), 63, 0.15, M);
                        TFD3 = Xckd(outputsigs(k,:), outputsigs(j,:), 1, 0.1, 0.1, M);
                        %% Plotting
                        %%% Thresholding
                        perc = 5;  % PLV Thresholding %
                        Q = length(SNR);
                        S0_m = cell(1,Q);
                        S1_m = cell(1,Q);
                        S2_m = cell(1,Q);
                        S3_m = cell(1,Q);
                        for q = 1:Q
                            S0_m{1,q} = S_stn{1,q}.*(abs(S_stn{1,q}) >= (100-perc*3)*max(max(abs(S_stn{1,q})))/100);
                            S1_m{1,q} = S_RID{1,q}.*(abs(S_RID{1,q}) >= (100-perc*3)*max(max(abs(S_RID{1,q})))/100);
                            S2_m{1,q} = S_WVD{1,q}.*(abs(S_WVD{1,q}) >= (100-perc)*max(max(abs(S_WVD{1,q})))/100);
                            S3_m{1,q} = S_CKD{1,q}.*(abs(S_CKD{1,q}) >= (100-perc)*max(max(abs(S_CKD{1,q})))/100);
                        end
                        %%% TFD Plotting
                        t = 0:1/fs:(N-1)/fs;
                        f = 0:fs/(2*M-1):fs/2;
                        tfr = 2;   % TFR plots resolution
                        figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                        subplot(2,2,1); flatwf(f,t(1:tfr:end),abs(TFD0(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(a) Ideal Time-Frequency Distribution','fontsize',10);
                        subplot(2,2,2); flatwf(f,t(1:tfr:end),abs(TFD2(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(b) Rihaczek Distribution','fontsize',10);
                        subplot(2,2,3); flatwf(f,t(1:tfr:end),abs(TFD1(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(c) Wigner-Ville Distribution','fontsize',10);
                        subplot(2,2,4); flatwf(f,t(1:tfr:end),abs(TFD3(:,1:tfr:end))','w','k');
                        xlabel('Frequency (Hz)'); ylabel('Time (s)'); title('(d) Compact Kernel Distribution (CKD)','fontsize',10);
                        %%% Brain Maps Plotting
                        for q = 1:length(SNR)
                            figure('Color',[1 1 1],'Position',[100, -100, 900 800]);
                            map = [(1:-0.01:0)',(1:-0.01:0)',(1:-0.01:0)']; colormap(map);
                            subplot(2,2,1), imagesc(abs(S0_m{1,q})); axis xy;
                            title('(a) Standard definition using Rihaczek','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,2), imagesc(abs(S1_m{1,q})); axis xy;
                            title('(b) Rihaczek Cross TFD','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,3), imagesc(abs(S2_m{1,q})); axis xy;
                            title('(c) Wigner-Ville Distribution','fontsize',10); xlabel('Channel'); ylabel('Channel');
                            subplot(2,2,4), imagesc(abs(S3_m{1,q})); axis xy;
                            title('(d) Compact Kernel Distribution (CKD)','fontsize',10); xlabel('Channel'); ylabel('Channel');
                        end
                        %% Saving Results
                        opt5 = input('Do you want to save the results (Y/N)\n','s');
                        if(opt5 == 'y' || opt5 == 'Y')
                            print(1,'PLV_TFD_2','-depsc');
                            print(2,'PLV_adult_sig2_50dB','-depsc');
                        else
                            return
                        end
                        ind2 = 0;
                    else
                        ind2 = 1;
                        fprintf(2,'Your selection must be either 1 or 2\n');
                    end
                end
                ind1 = 0;
            else
                ind1 = 1;
                fprintf(2,'Your selection must be either 1 or 2\n');
            end
        end
        ind = 0;
    else
        ind = 1;
        fprintf(2,'Your selection must be either 1 or 2\n');
    end
end