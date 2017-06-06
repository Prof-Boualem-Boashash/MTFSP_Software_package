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
fs        = 1;                          % Sampling frequency
N         = 256;                        % Number of samples
M         = 256;                        % Number of frequency bins
% Signal parameters
n         = 4;                          % Number of Source Signals
f_init    = [0.05, 0.2, 0.3, 0.45];     % LFM initial frequencies of n-1 sources
f_end     = [0.2, 0.05, 0.45, 0.3];     % LFM end frequencies of n-1 sources
% Reception parameters
m         = 3;                          % Number of Sensors
ra        = [15, 30, 45, 60];           % Reception angles of n sources
SNR       = 40;                         % Signal-to-noise Ratio
lambda    = 150;                        % Wavelength
d_space   = lambda/2;                   % Element spacing in the array of m antennas
% Reconstruction parameters
Eps_noise = 0.1;                        % remove points with energy less than Eps_noise/100
Eps_auto  = 0.9;                        % for removing cross-term
use_stft  = 1;                          % Use STFT in reconstruction

%% Computed variables
t  = 0:1/fs:(N-1)/fs;                   % Time array
f  = 0:fs/(2*M-1):fs/2;                 % Frequency array

%% Signal Generation
S = (signal_model(n, f_end, f_init, N, fs));
TFD_S = cell(1,n);
for i = 1:n
    TFD_S{1,i} = Xwvd(S(i,:), S(i,:), N-1, M);
end

%% Channel Mixing Model
rng(1);                   % Seed initialization of random noise
A = inst_model(n, m, ra, lambda, d_space); % Instantaneous Mixing Model
X = awgn(A*S, SNR);       % Additive White Gaussian Noise for provided SNR
TFD_X = cell(1,m);
for mm = 1:m
    TFD_X{1,mm} = Xspwvd(X(mm,:), X(mm,:), 'hann',31,'hann',31, M);
end

%% Whitening
W = whitening(X, n, m);

%% Multisensor Time-Frequency Distribution (MTFD)
D_wvd  = cell(m,m);
D_stft = cell(m,m);
D_mwvd = cell(m,m);
for m1 = 1:m
    for m2 = 1:m
        temp1 = Xwvd(X(m1,:), X(m2,:), N-1, M);
        D_wvd{m1,m2} = temp1.';
        if(m1 == m2)
            temp1 = stft(X(m1,:), M, hamming(51)');
            temp3 = (temp1'.^2).*D_wvd{m1,m2};
            D_stft{m1,m2} = temp1.';
            D_mwvd{m1,m2} = temp3;
        else
            D_stft{m1,m2} = zeros(N,M);
            D_mwvd{m1,m2} = zeros(N,M);
        end
    end
end

%% Auto-terms Selection
TFp = selatp(D_mwvd,'fast','auto',[Eps_noise, Eps_auto],W,D_wvd);

%% Concatenation of the WVD Matrix
TFM_xx = zeros(m*N, m*M);
for m1 = 1:m
    st1 = 1 + (m1-1)*M; fi1 = m1*M;
    for m2 = 1:m
        st2 = 1 + (m2-1)*N; fi2 = m2*N;
        TFM_xx(st2:fi2,st1:fi1) = D_wvd{m1,m2};
    end
end

%% STFD Eigen Vector Decomposition
Ntfp = size(TFp,1);
V  = zeros(Ntfp, 2*m);  % contains vectors for classification
Vw = zeros(Ntfp, 1);    % weights of vector in V
for ntfp = 1:Ntfp
    ti = TFp(ntfp,1); tt = ti + (0:m-1)*N;
    fi = TFp(ntfp,2); ff = fi + (0:m-1)*M;
    [Ub, Db, ~] = svd(TFM_xx(tt,ff));
    if(use_stft)
        temp1 = zeros(m,m);
        for m1 = 1:m
            temp1(m1,m1) = D_stft{m1,m1}(ti,fi);
        end
        temp1 = diag(temp1).*exp(-1j*angle(temp1(1,1)));
        temp2 = [real(temp1.') imag(temp1.')]/norm([real(temp1.') imag(temp1.')]);
        V(ntfp,:) = temp2;
    else
        u1 = Ub(:,1) * exp(-1j*angle(Ub(1,1)));
        V(ntfp,:) = [real(u1.') imag(u1.')];
    end
    Vw(ntfp) = Db(1,1)*Db(1,1)'; % eigenvalues
end
% sort the vectors in decreasing order of energy
Vinfo = [Vw V TFp];
Vinfo = sortrows(Vinfo,1); % increasing order
Vinfo = flipud(Vinfo);     % decreasing order
V     = Vinfo(:,2:2*m+1);
TFp   = Vinfo(:,2*m+2:2*m+3);

%% Clustering
map1 = kmeans(V, n,'Distance','sqEuclidean','Replicates',5,'EmptyAction','singleton');
Vmap = zeros(size(map1));
for k = 1:n
    ind = find(map1 == k);
    Vmap(ind,k) = k;
end

%% Sort TF signatures
win = 0;       % collect points in the window around a selected point
TFp_K = cell(1,n);
for k = 1:n
    ind = find(Vmap(:,k) == k);
    TFp_K{1,k} = TFp(ind,:);
end
TFD_r1 = getTFD(TFp_K, D_wvd, win);

%% Subspace Sepration of intersection points
[~, sig_order] = getIF(Vmap,TFD_r1,TFp);                % IF estimation
A_est = getmixture(Vmap,TFp,sig_order,D_stft,use_stft); % Mixture estimation
Vmap = zeros(size(Vmap));
TFD_r2 = cell(1,n);
for i = 1:n
    TFD_r2{1,i} = zeros(N,M);
end
Ntfp = size(TFp,1);
for i = 1:Ntfp
    ti = TFp(i,1); 
    fi = TFp(i,2); 
    C  = getSTFD(D_wvd,ti,fi);
    [U,DJ,~] = svd(C);
    V = U(:,1:2);
    P = (eye(m)-V*V');
    Q = P*A_est;
    for k = 1:size(Q,2)
        norm_Q(k) = norm(Q(:,k),'fro');
    end
    [~, idx] = sort(norm_Q,'ascend');
    A_ij = A_est(:,idx(1:2));
    B = pinv(A_ij);
    Dss = B*C*B';
    TFD_r2{1,idx(1)}(ti,fi) = real(Dss(1,1));
    TFD_r2{1,idx(2)}(ti,fi) = real(Dss(2,2));
end

%% Arrnaging Estimated signals
c1 = zeros(n,n);
c2 = zeros(n,n);
for i = 1:n
    for j = 1:n
        temp1 = TFD_r1{1,j};
        temp2 = TFD_r2{1,j};
        c1(i,j) = abs(corr(reshape(TFD_S{1,i},1,M*N)',reshape(temp1',1,M*N)'));
        c2(i,j) = abs(corr(reshape(TFD_S{1,i},1,M*N)',reshape(temp2',1,M*N)'));
    end
end
[~, I1] = max(c1);
[~, I2] = max(c2);
TFD_R1 = cell(1,n);
TFD_R2 = cell(1,n);
for i = 1:n
    TFD_R1{1,I1(i)} = TFD_r1{1,i}';
    TFD_R2{1,I2(i)} = TFD_r2{1,i}';
end

%% Plotting
tfr = 8;
figure('Color',[1 1 1],'Position',[100, -300, 1250 900]);
ha = tight_subplot(3,n,[0.05 0.01],[0.1 0.08],[0.07 0.07]);
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
    flatwf(f,t(1:tfr:end),abs(TFD_R1{1,i}(:,1:tfr:end))','w','k');
    title(['Estimated Source ' num2str(i)],'Fontsize',13);
    if i == 1,
        set(gca,'Xtick',[],'fontweight','bold','Fontsize',12);
        ylabel('Time (s)','Fontsize',12);
        legend('Clustering Method','location','SouthEast');
    elseif(i == 2)
        set(gca,'Ytick',[],'Xtick',[],'fontweight','bold','Fontsize',12);
        legend('Clustering Method','location','NorthEast');
    elseif(i == 3)
        set(gca,'Ytick',[],'Xtick',[],'fontweight','bold','Fontsize',12);
        legend('Clustering Method','location','NorthWest');
    elseif(i == 4)
        set(gca,'Ytick',[],'Xtick',[],'fontweight','bold','Fontsize',12);
        legend('Clustering Method','location','SouthWest');    
    end
    annotation('ellipse',[0.105673295454545 0.493975903614458 0.0344545454545455 0.0592921686746988],'Color',[1 0 0],'LineWidth',3);
    annotation('ellipse',[0.322036931818181 0.495316386188185 0.0344545454545455 0.0592921686746988],'Color',[1 0 0],'LineWidth',3);
    annotation('ellipse',[0.43112784090909 0.511402177072903 0.0344545454545455 0.0592921686746988],'Color',[1 0 0],'LineWidth',3);
    annotation('ellipse',[0.647491477272725 0.50469976420427 0.0344545454545453 0.0592921686746988],'Color',[1 0 0],'LineWidth',3);
    annotation('ellipse',[0.861127840909087 0.499337833909364 0.0344545454545453 0.0592921686746988],'Color',[1 0 0],'LineWidth',3);
    
    axes(ha(i+2*n));
    flatwf(f,t(1:tfr:end),abs(TFD_R2{1,i}(:,1:tfr:end))','w','k');
    title(['Estimated Source ' num2str(i)],'Fontsize',13);
    if i == 1,
        set(gca,'fontweight','bold','Fontsize',12);
        ylabel('Time (s)','Fontsize',12); xlabel('Frequency (Hz)','Fontsize',12);
        legend('Subspace Projection','location','SouthEast');
    elseif(i == 2)
        set(gca,'Ytick',[],'fontweight','bold','Fontsize',12);
        xlabel('Frequency (Hz)','Fontsize',12);
        legend('Subspace Projection','location','NorthEast');
    elseif(i == 3)
        set(gca,'Ytick',[],'fontweight','bold','Fontsize',12);
        xlabel('Frequency (Hz)','Fontsize',12);
        legend('Subspace Projection','location','NorthWest');
    elseif(i == 4)
        set(gca,'Ytick',[],'fontweight','bold','Fontsize',12);
        xlabel('Frequency (Hz)','Fontsize',12);
        legend('Subspace Projection','location','SouthWest');    
    end
end

tfr = 4;
figure('Color',[1 1 1],'Position',[100, 80, 500 450]);
ha = tight_subplot(1,1,[0.05 0.01],[0.12 0.08],[0.12 0.12]);
flatwf(f,t(1:tfr:end),abs(D_mwvd{1,1}(1:tfr:end,:)),'w','k');
ylabel('Time (s)','Fontsize',12); xlabel('Frequency (Hz)','Fontsize',12);
set(gca,'fontweight','bold');
title(['Received Signal on Sensor 1 (SNR = ' num2str(SNR) ' dB)'],'Fontsize',12);

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Inst_UBSS_Fig2','-depsc');
    print(2,'Inst_UBSS_Fig3','-depsc');
else
    return
end
