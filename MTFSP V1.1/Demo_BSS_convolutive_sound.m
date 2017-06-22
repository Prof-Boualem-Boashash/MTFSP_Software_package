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

clear; close all; clc;
addpath('Soundtracks');
addpath(genpath('Supporting functions'));

%% Loading signal
ind = 1;
while(ind)
    opt = input('Do you want to proceed with the default input signals (Y/N)\n','s');
    if(opt == 'y' || opt == 'Y')
        x1 = audioread('car_x.wav');
        x2 = audioread('seagull2.wav');
        ind = 0;
    elseif(opt == 'n' || opt == 'N')
        ind1 = 1; ind2 = 1;
        disp('You can select a signal by entering its corrosponding number');
        disp('-------------------------------------------------------------');
        fprintf('1. Bird Song\n2. Bird Chirp\n3. Car Startup\n4. Dolphine Click\n5. Jet Aircraft\n6. Seagulls\n');
        disp('-------------------------------------------------------------');
        while(ind1)
            opt1 = input('Kindly select the first signal  : ');
            switch opt1
                case 1,  ind1 = 0; x1 = audioread('bird.wav');
                case 2,  ind1 = 0; x1 = audioread('bird_chirp.wav');
                case 3,  ind1 = 0; x1 = audioread('car_x.wav');
                case 4,  ind1 = 0; x1 = audioread('dolphin.wav');
                case 5,  ind1 = 0; x1 = audioread('jet_doppler2.wav');
                case 6,  ind1 = 0; x1 = audioread('seagull2.wav');
                otherwise, ind1 = 1; fprintf(2,'Your selection must be betwen 1 and 6\n');
            end
        end
        while(ind2)
            opt2 = input('Kindly select the second signal : ');
            switch opt2
                case 1,  ind2 = 0; x2 = audioread('bird.wav');
                case 2,  ind2 = 0; x2 = audioread('bird_chirp.wav');
                case 3,  ind2 = 0; x2 = audioread('car_x.wav');
                case 4,  ind2 = 0; x2 = audioread('dolphin.wav');
                case 5,  ind2 = 0; x2 = audioread('jet_doppler2.wav');
                case 6,  ind2 = 0; x2 = audioread('seagull2.wav');
                otherwise, ind2 = 1; fprintf(2,'Your selection must be betwen 1 and 6\n');
            end
        end
        ind = 0;
    else
        ind = 1;
        fprintf(2,'Your selection must be either [Y] or [N]\n');
    end
end

%% Computed Variables
fs = 11025;             % Sampling Frequency
N  = length(x1);        % Signal length
t  = 0:1/fs:(N-1)/fs;   % Time array
f  = 0:fs/(2*N-1):fs/2; % Frequency Array
Q  = 1024;              % Number of STFD samples
M  = 1024;              % Number of frequency bins
n  = 2;                 % Number of Source Signals
m  = 3;                 % Number of Sensors

%% Signal Generation
x1 = x1./max(x1); x2 = x2./max(x2);
S = [x1 x2]';
TFD_S = cell(1,n);
for i = 1:n
    TFD_S{1,i} = abs(spectrogram(S(i,:),hamming(127),63,1024,fs));
end
[~,ffs,tts] = spectrogram(S(1,:),hamming(127),63,1024,fs);

%% Channel Mixing Model
L = 1; LL = 2;
syms z1
Az = [1.0 + 0.0*z1  0.85 + 0.10*z1
      0.7 + 0.4*z1  0.25 + 1.00*z1
      1.0 + 0.5*z1  0.70 + 0.85*z1];
[X, A] = conv_model(S, Az, L, LL);
TFD_X = cell(1,n*(LL+L));
for i = 1:n*(LL+L)
    TFD_X{1,i} = abs(spectrogram(X(i,:),hamming(127),63,1024,fs));
end

%% Whitening
W = whitening(X, n*(LL+L), m*LL);
Z = W*X;

%% Source Separation
% Multisensor Time-Frequency Distribution (MTFD)
QQ = 1000;
zz = Z(:,QQ:Q+QQ-1);
D = mtfd(zz, 'wvd', Q-LL-1, M);
% Selection of Auto-STFDs
Ns = 6;                    % Number of Auto-terms matrices used in JD
th = 1e-1;                 % Thereshold for Auto-terms and Cross-terms selection
aSTFD_all = select_TFD_Convolutive(D, n, LL+L, th);
aSTFD =  zeros(n*(LL+L), n*(LL+L),Ns);
Mx = floor(size(aSTFD_all,3)/Ns);
for ii = 1:Ns
    aSTFD(:,:,ii) = real(mean(aSTFD_all(:,:,(ii-1)*Mx+1:ii*Mx),3));
end
% Joint Block Diagonalization
V = JointBlockDiag(aSTFD, (LL+L)*ones(1,n),1e-9,5000);
R = V'*Z;

%% Arrnaging Estimated signals
TFD_R = cell(1,n*(LL+L));
for i = 1:n*(LL+L)
    R(i,:) = R(i,:)./max(abs(R(i,:)));
    TFD_R{1,i} = abs(spectrogram(R(i,:),hamming(127),63,1024,fs));
end
c = zeros(n,n*(LL+L));
for i = 1:n
    for j = 1:n*(LL+L)
        [Ms, Ns] = size(TFD_S{1,i}); [Mr, Nr] = size(TFD_R{1,j});
        c(i,j) = abs(corr(reshape(TFD_S{1,i},1,Ms*Ns)',reshape(TFD_R{1,j},1,Mr*Nr)'));
    end
end
[~, I] = sort(c,2,'descend');
I = I(:,1:(LL+L));
I = reshape(I,1,6);
R = R(I,:);
for i = 1:n*(LL+L)
    TFD_R{1,i} = abs(spectrogram(R(i,:),hamming(127),63,1024,fs));
end
[~,ffr,ttr] = spectrogram(R(1,:),hamming(127),63,1024,fs);

%% Plotting
figure('Color',[1 1 1],'Position',[100, 80, 750 600]); cnt = 0;
ha = tight_subplot(n*(LL+L)/2+1,n,[0.05 0.01],[0.1 0.07],[0.08 0.08]);
for i = 1:n*(LL+L+1)
    if i <= n
        axes(ha(i));
        imagesc(ffs,tts,log10(TFD_S{1,i}')); axis xy; title(['Source Signal ' num2str(i)],'Fontsize',12);
        if i == 1,
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Time (s)','Fontsize',12);
        else
            set(gca,'Ytick',[],'Xtick',[]);
        end
    else
        axes(ha(i));
        imagesc(ffr,ttr,log10(TFD_R{1,i-n}')); axis xy;
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

figure('Color',[1 1 1],'Position',[100, 80, 750 600]); cnt = 0;
ha = tight_subplot(n*(LL+L)/2+1,n,[0.05 0.01],[0.1 0.07],[0.08 0.08]);
for i = 1:n*(LL+L+1)
    if i <= n
        axes(ha(i));
        plot(t,S(i,:)); title(['Source Signal ' num2str(i)],'Fontsize',12); ylim([-1 1]);
        if i == 1,
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Amplitude','Fontsize',12);
        else
            set(gca,'Ytick',[],'Xtick',[]);
        end
    else
        axes(ha(i));
        plot(t(1:end-LL),R(i-n,:)); ylim([-1 1])
        if mod(i,2) && i ~= n*(LL+L+1)-1
            cnt = cnt + 1;
            set(gca,'Xtick',[],'fontweight','bold');
            ylabel('Amplitude','Fontsize',12);
            title(['Estimated Source 1 v.' num2str(cnt)],'Fontsize',12)
        elseif i == n*(LL+L+1)-1
            cnt = cnt + 1;
            set(gca,'fontweight','bold');
            ylabel('Amplitude','Fontsize',12); xlabel('Time (s)','Fontsize',12);
            title(['Estimated Source 1 v.' num2str(cnt)],'Fontsize',12)
        elseif i == n*(LL+L+1)
            set(gca,'Ytick',[],'fontweight','bold');
            xlabel('Time (s)','Fontsize',12);
            title(['Estimated Source 2 v.' num2str(cnt)],'Fontsize',12)
        else
            set(gca,'Ytick',[],'Xtick',[]);
            title(['Estimated Source 2 v.' num2str(cnt)],'Fontsize',12)
        end
    end
end

figure('Color',[1 1 1],'Position',[100, 80, 500 450]);
ha = tight_subplot(1,1,[0.05 0.01],[0.12 0.08],[0.12 0.12]);
imagesc(ffr,ttr,log10(TFD_X{1,4}')); axis xy;
ylabel('Time (s)','Fontsize',12); xlabel('Frequency (Hz)','Fontsize',12);
set(gca,'fontweight','bold'); title('Received Signal on Sensor 4','Fontsize',14);

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'Conv_BSS_Sound_Fig1','-depsc');
    print(2,'Conv_BSS_Sound_Fig2','-depsc');
    print(3,'Conv_BSS_Sound_Fig3','-depsc');
else
    return
end