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
        %% Extracting SNR information from Data Folders
        S = dir('DOA Data');
        SNR_N = sum([S(~ismember({S.name},{'.','..'})).isdir]);
        SNR_i = zeros(1,SNR_N);
        for i = 1:SNR_N
            [~, r1] = strtok(S(7+i).name);
            [s1, r2] = strtok(r1);
            SNR_i(1,i) = str2double(s1);
        end
        SNR_i = sort(SNR_i);
        
        %% MUSIC and TF-MUSIC Spectrum
        SNR_N = 2;
        figure('Color',[1 1 1],'Position',[100, 70, 900 500]);
        ha = tight_subplot(SNR_N/2,2,[0.05 0.01],[0.12 0.12],[0.08 0.08]);
        for i = 1:SNR_N
            Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
            load([Path '\Averaged_Spectrum']);
            axes(ha(i));
            plot(theta, P_tf_music_avg,'-b','linewidth',2); hold on;
            plot(theta, P_music_avg,'-.r','linewidth',2); hold on;
            plot(repmat(ra(1),10),linspace(0,1.2,10),'k--','linewidth',2); hold on
            plot(repmat(ra(2),10),linspace(0,1.2,10),'k--','linewidth',2); grid on
            axis([0 39.9 0 1.2])
            if(mod(i,2) && i < SNR_N-1)
                set(gca,'XTickLabel','','fontweight','bold','fontsize',13);
                ylabel('P_M_U_S_I_C (\theta)','Fontsize',16);
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            elseif i == SNR_N-1
                set(gca,'fontweight','bold','fontsize',13);
                ylabel('P_M_U_S_I_C (\theta)','Fontsize',16);
                xlabel('\theta (deg)','Fontsize',16);
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            elseif i == SNR_N
                set(gca,'YTickLabel','','fontweight','bold','fontsize',13);
                xlabel('\theta (deg)','Fontsize',16);
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
                legend('TF-MUSIC Averaged Spectrum','MUSIC Averaged Spectrum','True Angles','Location','Southwest')
            else
                set(gca,'XTickLabel','','YTickLabel','','fontweight','bold');
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            end
        end
        set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
        
        %% DOA NMSE
        SNR_N = sum([S(~ismember({S.name},{'.','..'})).isdir]);
        Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
        load([Path '\DOA'],'ra');
        nmse_music     = zeros(length(ra), SNR_N);
        nmse_tf_music  = zeros(length(ra), SNR_N);
        nmse_esprit    = zeros(length(ra), SNR_N);
        nmse_tf_esprit = zeros(length(ra), SNR_N);
        music_rate    = zeros(length(ra), SNR_N);
        tf_music_rate = zeros(length(ra), SNR_N);
        for i = 1:SNR_N
            Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
            load([Path '\DOA']);
            for j = 1:length(ra)
                temp1 = DOA_music(:,j); temp2 = DOA_tf_music(:,j);
                music_rate(j,i)    = sum(temp1~=0)/length(temp1);
                tf_music_rate(j,i) = sum(temp2~=0)/length(temp2);
                temp1 = temp1(temp1~=0); temp2 = temp2(temp2~=0);
                nmse_music(j,i)     = (mean(((temp1 - ra(j))/ra(j)).^2));
                nmse_tf_music(j,i)  = (mean(((temp2 - ra(j))/ra(j)).^2));
                nmse_esprit(j,i)    = (mean(((DOA_esprit(:,j) - ra(j))/ra(j)).^2));
                nmse_tf_esprit(j,i) = (mean(((DOA_tf_esprit(:,j) - ra(j))/ra(j)).^2));
            end
        end
        fprintf(2,'Mean Probability of Detection (Pd)\n');
        fprintf('MUSIC     : %0.3f\n',mean(mean(music_rate)));
        fprintf('TF_MUSIC  : %0.3f\n',mean(mean(tf_music_rate)))
        fprintf('ESPRIT    : %0.1f\n',1)
        fprintf('TF_ESPRIT : %0.1f\n',1)
        figure('Color',[1 1 1],'Position',[100, 10, 650, 550]);
        ha = tight_subplot(1,1,[0.01 0.01],[0.12 0.1],[0.12 0.12]);
        axes(ha(1));
        plot(SNR_i,10*log10(mean(nmse_tf_music)),'-b','linewidth',2);  hold on;
        plot(SNR_i,10*log10(mean(nmse_music)),'--b','linewidth',2);    hold on;
        plot(SNR_i,10*log10(mean(nmse_tf_esprit)),'-.r','linewidth',2); hold on;
        plot(SNR_i,10*log10(mean(nmse_esprit)),':r','linewidth',2);   grid on;
        xlim([SNR_i(1) SNR_i(end)]);
        legend('TF-MUSIC','MUSIC','TF-ESPRIT','ESPRIT','Location','Southwest');
        set(gca,'fontweight','bold','fontsize',14);
        title('DOA Normalized Mean Square Error','fontsize',18);
        xlabel('SNR (dB)'); ylabel('NMSE (dB)');
        set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
        
        %% DOA PDFs
        SNR_N = 2;
        bins_N = 100;
        for i = 1:SNR_N
            aa= figure('Color',[1 1 1],'Position',[100, 0, 650, 700]);
            ha = tight_subplot(2,1,[0.01 0.01],[0.1 0.1],[0.05 0.05]);
            Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
            load([Path '\DOA']);
            axes(ha(1)); temp = DOA_music; temp = temp(temp~=0);
            histogram(temp, bins_N,'Normalization','probability','facecolor',[1 0.84 0],'linewidth',1.2); hold on;
            temp = DOA_tf_music; temp = temp(temp~=0);
            histogram(temp, bins_N,'Normalization','probability',...
                'facecolor',[1 0 0],'linewidth',1.2); xlim([0 40]);
            temp1 = DOA_music(:,1); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_music(:,2); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag1 = ['         MUSIC' char(10) tag11 char(10) tag12];
            temp1 = DOA_tf_music(:,1); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_tf_music(:,2); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag2 = ['       TF-MUSIC' char(10) tag11 char(10) tag12];
            legend(tag1, tag2,'Location','Northwest'); grid on;
            set(gca,'XTickLabel','','YTickLabel','','fontweight','bold','fontsize',13);
            title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            axes(ha(2));
            histogram(DOA_esprit, bins_N,'Normalization','probability',...
                'facecolor',[0 0.5 0],'linewidth',1.2); hold on;
            histogram(DOA_tf_esprit, bins_N,'Normalization','probability',...
                'facecolor',[0 0.45 0.74],'linewidth',1.2); xlim([0 40]);
            temp1 = DOA_esprit(:,1); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_esprit(:,2); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag1 = ['         ESPRIT' char(10) tag11 char(10) tag12];
            temp1 = DOA_tf_esprit(:,1); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_tf_esprit(:,2); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag2 = ['       TF-ESPRIT' char(10) tag11 char(10) tag12];
            legend(tag1, tag2,'Location','Northwest'); grid on;
            set(gca,'YTickLabel','','fontweight','bold','fontsize',13);
            xlabel('\theta (deg)');
            set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
            set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
        end
        
        %% Saving
        opt1 = input('Do you want to save the results (Y/N)\n','s');
        if(opt1 == 'y' || opt1 == 'Y')
            print(1,'DOA_Spectrum','-dpdf','-r1024');
            print(2,'DOA_NMSE','-dpdf','-r1024');
            for i = 3:SNR_N+2
                print(i,['PDF_SNR(' num2str(SNR_i(i-2)),')'],'-dpdf','-r1024');
            end
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
            disp(['SNR ' num2str(SNR_i(i)), ' dB'])
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
            Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB']; mkdir(Path);
            save([Path '/Averaged_Spectrum'],'P_music_avg','P_tf_music_avg','theta','ra','-v7.3');
            save([Path '/DOA'],'DOA_music','DOA_tf_music','DOA_esprit','DOA_tf_esprit','SNR','ra','-v7.3')
        end
        disp('Finished')
        %% Extracting SNR information from Data Folders
        S = dir('DOA Data');
        SNR_N = sum([S(~ismember({S.name},{'.','..'})).isdir]);
        SNR_i = zeros(1,SNR_N);
        for i = 1:SNR_N
            [~, r1] = strtok(S(5+i).name);
            [~, r2] = strtok(r1);
            [s1, r3] = strtok(r2);
            SNR_i(1,i) = str2double(s1);
        end
        SNR_i = sort(SNR_i);
        
        %% MUSIC and TF-MUSIC Spectrum
        SNR_N = 2;
        figure('Color',[1 1 1],'Position',[100, 70, 900 500]);
        ha = tight_subplot(SNR_N/2,2,[0.05 0.01],[0.12 0.12],[0.08 0.08]);
        for i = 1:SNR_N
            Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
            load([Path '\Averaged_Spectrum']);
            axes(ha(i));
            plot(theta, P_tf_music_avg,'-b','linewidth',2); hold on;
            plot(theta, P_music_avg,'-.r','linewidth',2); hold on;
            plot(repmat(ra(1),10),linspace(0,1.2,10),'k--','linewidth',2); hold on
            plot(repmat(ra(2),10),linspace(0,1.2,10),'k--','linewidth',2); grid on
            axis([0 39.9 0 1.2])
            if(mod(i,2) && i < SNR_N-1)
                set(gca,'XTickLabel','','fontweight','bold','fontsize',13);
                ylabel('P_M_U_S_I_C (\theta)','Fontsize',16);
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            elseif i == SNR_N-1
                set(gca,'fontweight','bold','fontsize',13);
                ylabel('P_M_U_S_I_C (\theta)','Fontsize',16);
                xlabel('\theta (deg)','Fontsize',16);
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            elseif i == SNR_N
                set(gca,'YTickLabel','','fontweight','bold','fontsize',13);
                xlabel('\theta (deg)','Fontsize',16);
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
                legend('TF-MUSIC Averaged Spectrum','MUSIC Averaged Spectrum','True Angles','Location','Southwest')
            else
                set(gca,'XTickLabel','','YTickLabel','','fontweight','bold');
                title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            end
        end
        set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
        
        %% DOA NMSE
        SNR_N = sum([S(~ismember({S.name},{'.','..'})).isdir]);
        Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
        load([Path '\DOA'],'ra');
        nmse_music     = zeros(length(ra), SNR_N);
        nmse_tf_music  = zeros(length(ra), SNR_N);
        nmse_esprit    = zeros(length(ra), SNR_N);
        nmse_tf_esprit = zeros(length(ra), SNR_N);
        music_rate    = zeros(length(ra), SNR_N);
        tf_music_rate = zeros(length(ra), SNR_N);
        for i = 1:SNR_N
            Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
            load([Path '\DOA']);
            for j = 1:length(ra)
                temp1 = DOA_music(:,j); temp2 = DOA_tf_music(:,j);
                music_rate(j,i)    = sum(temp1~=0)/length(temp1);
                tf_music_rate(j,i) = sum(temp2~=0)/length(temp2);
                temp1 = temp1(temp1~=0); temp2 = temp2(temp2~=0);
                nmse_music(j,i)     = (mean(((temp1 - ra(j))/ra(j)).^2));
                nmse_tf_music(j,i)  = (mean(((temp2 - ra(j))/ra(j)).^2));
                nmse_esprit(j,i)    = (mean(((DOA_esprit(:,j) - ra(j))/ra(j)).^2));
                nmse_tf_esprit(j,i) = (mean(((DOA_tf_esprit(:,j) - ra(j))/ra(j)).^2));
            end
        end
        fprintf(2,'Mean Probability of Detection (Pd)\n');
        fprintf('MUSIC     : %0.3f\n',mean(mean(music_rate)));
        fprintf('TF_MUSIC  : %0.3f\n',mean(mean(tf_music_rate)))
        fprintf('ESPRIT    : %0.1f\n',1)
        fprintf('TF_ESPRIT : %0.1f\n',1)
        figure('Color',[1 1 1],'Position',[100, 10, 650, 550]);
        ha = tight_subplot(1,1,[0.01 0.01],[0.12 0.1],[0.12 0.12]);
        axes(ha(1));
        plot(SNR_i,10*log10(mean(nmse_tf_music)),'-b','linewidth',2);  hold on;
        plot(SNR_i,10*log10(mean(nmse_music)),'--b','linewidth',2);    hold on;
        plot(SNR_i,10*log10(mean(nmse_tf_esprit)),'-.r','linewidth',2); hold on;
        plot(SNR_i,10*log10(mean(nmse_esprit)),':r','linewidth',2);   grid on;
        xlim([SNR_i(1) SNR_i(end)]);
        legend('TF-MUSIC','MUSIC','TF-ESPRIT','ESPRIT','Location','Southwest');
        set(gca,'fontweight','bold','fontsize',14);
        title('DOA Normalized Mean Square Error','fontsize',18);
        xlabel('SNR (dB)'); ylabel('NMSE (dB)');
        set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
        set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
        
        %% DOA PDFs
        SNR_N = 2;
        bins_N = 100;
        for i = 1:SNR_N
            aa= figure('Color',[1 1 1],'Position',[100, 0, 650, 700]);
            ha = tight_subplot(2,1,[0.01 0.01],[0.1 0.1],[0.05 0.05]);
            Path = ['DOA Data\SNR ' num2str(SNR_i(i)), ' dB'];
            load([Path '\DOA']);
            axes(ha(1)); temp = DOA_music; temp = temp(temp~=0);
            histogram(temp, bins_N,'Normalization','probability','facecolor',[1 0.84 0],'linewidth',1.2); hold on;
            temp = DOA_tf_music; temp = temp(temp~=0);
            histogram(temp, bins_N,'Normalization','probability',...
                'facecolor',[1 0 0],'linewidth',1.2); xlim([0 40]);
            temp1 = DOA_music(:,1); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_music(:,2); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag1 = ['         MUSIC' char(10) tag11 char(10) tag12];
            temp1 = DOA_tf_music(:,1); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_tf_music(:,2); temp1 = temp1(temp1~=0); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag2 = ['       TF-MUSIC' char(10) tag11 char(10) tag12];
            legend(tag1, tag2,'Location','Northwest'); grid on;
            set(gca,'XTickLabel','','YTickLabel','','fontweight','bold','fontsize',13);
            title(['SNR = ' num2str(SNR_i(i)), ' dB'],'fontsize',20);
            axes(ha(2));
            histogram(DOA_esprit, bins_N,'Normalization','probability',...
                'facecolor',[0 0.5 0],'linewidth',1.2); hold on;
            histogram(DOA_tf_esprit, bins_N,'Normalization','probability',...
                'facecolor',[0 0.45 0.74],'linewidth',1.2); xlim([0 40]);
            temp1 = DOA_esprit(:,1); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_esprit(:,2); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag1 = ['         ESPRIT' char(10) tag11 char(10) tag12];
            temp1 = DOA_tf_esprit(:,1); u1 = mean(temp1); s = std(temp1);
            tag11 = ['\mu_1 = ', num2str(round(u1,1)), ', \sigma_1 = ', num2str(round(s,1))];
            temp1 = DOA_tf_esprit(:,2); u1 = mean(temp1); s = std(temp1);
            tag12 = ['\mu_2 = ', num2str(round(u1,1)), ', \sigma_2 = ', num2str(round(s,1))];
            tag2 = ['       TF-ESPRIT' char(10) tag11 char(10) tag12];
            legend(tag1, tag2,'Location','Northwest'); grid on;
            set(gca,'YTickLabel','','fontweight','bold','fontsize',13);
            xlabel('\theta (deg)');
            set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
            set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
        end
        
        %% Saving
        opt1 = input('Do you want to save the results (Y/N)\n','s');
        if(opt1 == 'y' || opt1 == 'Y')
            print(1,'DOA_Spectrum','-dpdf','-r1024');
            print(2,'DOA_NMSE','-dpdf','-r1024');
            for i = 3:SNR_N+2
                print(i,['PDF_SNR(' num2str(SNR_i(i-2)),')'],'-dpdf','-r1024');
            end
        else
            return
        end
        ind = 0;
    else
        ind = 1;
        fprintf(2,'Your selection must be either 1 or 2\n');
    end
end