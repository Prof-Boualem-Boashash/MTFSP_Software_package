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
addpath(genpath('Lead Field Matrix Data'));
addpath(genpath('DOA Data'));

%% Main
ind = 1;
while(ind)
    fprintf('Please select one EEG from the following:\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('1. Newborn ');fprintf(2,'(Only Presenting results)\n');
    fprintf('2. Adult\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('You can select an option by entering its corrosponding number 1 or 2\n');
    opt = input('');
    if(opt == 1 && isnumeric(opt))
        %% Loading Results
        run Combine_Data.m;
        load('LeadfieldMatrix_adult.mat','Brain');
        load DOA_EEG_newborn.mat
        n = 3;
        orientation = [0.5 0.9 0.8; -0.026 -0.36 -0.46 ; -0.24 -0.16 -0.2];
        
        %% Plotting
        alpha = 0.4; G = 5;
        figure('Color',[1 1 1],'Position',[100, 0, 900 600]);
        ha = tight_subplot(2,n,[0.05 0.01],[0.01 0.1],[0.01 0.01]);
        for i = 1:n
            axes(ha(i));
            scatter3(Brain(:,1),Brain(:,2),Brain(:,3),'k.','MarkerEdgeAlpha',alpha,'sizedata',G); hold on
            plot3(Brain(DOA_music,1),Brain(DOA_music,2),Brain(DOA_music,3),'r.','markersize',15); hold on
            plot3(Brain(locations,1),Brain(locations,2),Brain(locations,3),'ko','markersize',10,'linewidth',3);
            axis off; set(gca,'fontsize',10,'fontweight','bold');
            if(i==1), view(0,90); legend('Top View','location','NorthEast');
            elseif(i==2), view(90,0); legend('Back View','location','NorthEast');
                title('Conventional MUSIC','fontsize',16);
            else view(0,0); legend('Side View','location','NorthEast');
            end
        end
        for i = 1+n:2*n
            axes(ha(i));
            scatter3(Brain(:,1),Brain(:,2),Brain(:,3),'k.','MarkerEdgeAlpha',alpha,'sizedata',G); hold on
            plot3(Brain(DOA_tf_music,1),Brain(DOA_tf_music,2),Brain(DOA_tf_music,3),'r.','markersize',15); hold on
            plot3(Brain(locations,1),Brain(locations,2),Brain(locations,3),'ko','markersize',10,'linewidth',3);
            axis off; set(gca,'fontsize',10,'fontweight','bold');
            if(i==1+n), view(0,90); legend('Top View','location','NorthEast');
            elseif(i==2+n), view(90,0); legend('Back View','location','NorthEast');
                title('Time-Frequency MUSIC','fontsize',16);
            else view(0,0); legend('Side View','location','NorthEast');
            end
        end
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'PaperOrientation','landscape');
        
        %% Saving
        opt1 = input('Do you want to save the results (Y/N)\n','s');
        if(opt1 == 'y' || opt1 == 'Y')
            print('DOA_EEG_newborn','-dpdf','-r512');
        else
            return
        end
        ind = 0;
        
    elseif(opt == 2 && isnumeric(opt))
        ind1 = 1;
        while(ind1)
            fprintf('Please select one simulation option from the following:\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('1. Present pre-computed results\n2. Re-compute and present results\n');
            fprintf('-------------------------------------------------------------\n');
            fprintf('You can select an option by entering its corrosponding number 1 or 2\n');
            opt1 = input('');
            if(opt1 == 1 && isnumeric(opt1))
                %% Loading Results
                run Combine_Data.m;
                load('LeadfieldMatrix_adult.mat','Brain');
                load DOA_EEG_adult.mat
                n = 3;
                orientation = [0.5 0.9 0.8; -0.026 -0.36 -0.46 ; -0.24 -0.16 -0.2];
                
                %% Plotting
                alpha = 0.4; G = 5;
                figure('Color',[1 1 1],'Position',[100, 0, 900 600]);
                ha = tight_subplot(2,n,[0.05 0.01],[0.01 0.1],[0.01 0.01]);
                for i = 1:n
                    axes(ha(i));
                    scatter3(Brain(:,1),Brain(:,2),Brain(:,3),'k.','MarkerEdgeAlpha',alpha,'sizedata',G); hold on
                    plot3(Brain(DOA_music,1),Brain(DOA_music,2),Brain(DOA_music,3),'r.','markersize',15); hold on
                    plot3(Brain(locations,1),Brain(locations,2),Brain(locations,3),'ko','markersize',10,'linewidth',3);
                    axis off; set(gca,'fontsize',10,'fontweight','bold');
                    if(i==1), view(0,90); legend('Top View','location','NorthEast');
                    elseif(i==2), view(90,0); legend('Back View','location','NorthEast');
                        title('Conventional MUSIC','fontsize',16);
                    else view(0,0); legend('Side View','location','NorthEast');
                    end
                end
                for i = 1+n:2*n
                    axes(ha(i));
                    scatter3(Brain(:,1),Brain(:,2),Brain(:,3),'k.','MarkerEdgeAlpha',alpha,'sizedata',G); hold on
                    plot3(Brain(DOA_tf_music,1),Brain(DOA_tf_music,2),Brain(DOA_tf_music,3),'r.','markersize',15); hold on
                    plot3(Brain(locations,1),Brain(locations,2),Brain(locations,3),'ko','markersize',10,'linewidth',3);
                    axis off; set(gca,'fontsize',10,'fontweight','bold');
                    if(i==1+n), view(0,90); legend('Top View','location','NorthEast');
                    elseif(i==2+n), view(90,0); legend('Back View','location','NorthEast');
                        title('Time-Frequency MUSIC','fontsize',16);
                    else view(0,0); legend('Side View','location','NorthEast');
                    end
                end
                set(gcf,'PaperPositionMode','auto');
                set(gcf,'PaperOrientation','landscape');
                
                %% Saving
                opt2 = input('Do you want to save the results (Y/N)\n','s');
                if(opt2 == 'y' || opt2 == 'Y')
                    print('DOA_EEG_adult','-dpdf','-r512');
                else
                    return
                end
                ind1 = 0;
                
            elseif(opt1 == 2 && isnumeric(opt1))
                %% Parameters that can be changed by the user
                fs   = 20;   % Sampling Frequency
                T    = 12;   % Total Duration
                SNR  = 3.5;  % Signal-to-Noise Level in dB
                perc = 0.1;  % Percentage of the STFD maximum power to select high-energy (t,f) points
                K    = 1000; % Number of Iterations
                rng(3);
                
                %% Signal Generation
                t = 1:1/fs:(T*fs-1)/fs;
                sigs(1,:) = gauspuls(t-T/5,5,0.2);
                sigs(2,:) = gauspuls(t-T/2,6,0.3);
                tmp = chirp(t(1:ceil(length(t)*0.1)),1,4,5);
                sigs(3,ceil(length(t)*0.65):ceil(length(t)*0.65)-1+length(tmp)) = tmp;
                n = size(sigs,1);
                N = length(sigs);
                M = 2^nextpow2(N);
                
                %% main
                DOA_music = zeros(K, n);
                DOA_tf_music = zeros(K, n);
                run Combine_Data.m;
                load LeadfieldMatrix_adult.mat;
                ch_n   = size(Gain,1);
                [A, B] = size(Brain);
                orientation = [0.5 0.9 0.8; -0.026 -0.36 -0.46 ; -0.24 -0.16 -0.2];
                locations = [32836 300 15050];
                for k = 1:K
                    outputsigs = zeros(ch_n, N);
                    for i = 1:n
                        gains = Gain(:,n*locations+i-1).*(repmat(orientation(i,:),ch_n,1));
                        outputsigs = outputsigs + gains*sigs;
                    end
                    noisyOut = zeros(ch_n, N);
                    for i = 1:ch_n
                        noisyOut(i,:) = awgn(outputsigs(i,:),SNR,'measured');
                    end
                    %% Time-Domain MUSIC
                    DOA_music(k,:) = EEG_music(noisyOut, n, ch_n, A, Gain, orientation);
                    %% Time-Frequency MUSIC
                    D = mtfd(noisyOut, 'wvd', N/4, M);
                    %Averaged Auto-TFD
                    D_avg = zeros(M, N);
                    for mm = 1:ch_n, D_avg = D{mm,mm} + D_avg; end
                    D_avg = D_avg./ch_n;
                    %Selection of high-energy (t,f) points
                    thr = perc*max(max(D_avg));
                    Tr = abs(D_avg) >= thr;
                    [F_trace, ~] = find(Tr);
                    n_p = length(F_trace);
                    D_s = zeros(ch_n, ch_n);
                    for m1 = 1:ch_n
                        for m2 = 1:ch_n
                            D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
                        end
                    end
                    clear D;
                    DOA_tf_music(k,:) = EEG_tf_music(D_s, n, ch_n, A, Gain, orientation);
                    disp(100*k/K)
                end
                clear Gain
                
                %% Saving
                save('DOA Data\DOA_EEG_adult','DOA_music','DOA_tf_music','SNR','locations','-v7.3');
                
                %% Loading Results
                run Combine_Data.m;
                load('LeadfieldMatrix_adult.mat','Brain');
                load DOA_EEG_adult.mat
                n = 3;
                orientation = [0.5 0.9 0.8; -0.026 -0.36 -0.46 ; -0.24 -0.16 -0.2];
                
                %% Plotting
                alpha = 0.4; G = 5;
                figure('Color',[1 1 1],'Position',[100, 0, 900 600]);
                ha = tight_subplot(2,n,[0.05 0.01],[0.01 0.1],[0.01 0.01]);
                for i = 1:n
                    axes(ha(i));
                    scatter3(Brain(:,1),Brain(:,2),Brain(:,3),'k.','MarkerEdgeAlpha',alpha,'sizedata',G); hold on
                    plot3(Brain(DOA_music,1),Brain(DOA_music,2),Brain(DOA_music,3),'r.','markersize',15); hold on
                    plot3(Brain(locations,1),Brain(locations,2),Brain(locations,3),'ko','markersize',10,'linewidth',3);
                    axis off; set(gca,'fontsize',10,'fontweight','bold');
                    if(i==1), view(0,90); legend('Top View','location','NorthEast');
                    elseif(i==2), view(90,0); legend('Back View','location','NorthEast');
                        title('Conventional MUSIC','fontsize',16);
                    else view(0,0); legend('Side View','location','NorthEast');
                    end
                end
                for i = 1+n:2*n
                    axes(ha(i));
                    scatter3(Brain(:,1),Brain(:,2),Brain(:,3),'k.','MarkerEdgeAlpha',alpha,'sizedata',G); hold on
                    plot3(Brain(DOA_tf_music,1),Brain(DOA_tf_music,2),Brain(DOA_tf_music,3),'r.','markersize',15); hold on
                    plot3(Brain(locations,1),Brain(locations,2),Brain(locations,3),'ko','markersize',10,'linewidth',3);
                    axis off; set(gca,'fontsize',10,'fontweight','bold');
                    if(i==1+n), view(0,90); legend('Top View','location','NorthEast');
                    elseif(i==2+n), view(90,0); legend('Back View','location','NorthEast');
                        title('Time-Frequency MUSIC','fontsize',16);
                    else view(0,0); legend('Side View','location','NorthEast');
                    end
                end
                set(gcf,'PaperPositionMode','auto');
                set(gcf,'PaperOrientation','landscape');
                
                %% Saving
                opt2 = input('Do you want to save the results (Y/N)\n','s');
                if(opt2 == 'y' || opt2 == 'Y')
                    print('DOA_EEG_adult','-dpdf','-r512');
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
    else
        ind = 1;
        fprintf(2,'Your selection must be either 1 or 2\n');
    end
end
