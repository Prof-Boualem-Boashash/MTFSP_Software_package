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

%% Parameters that can be changed by the user
% General parameters
fs      = 1;                     % Sampling frequency
N       = 512;                   % Number of samples
M       = 512;                   % Number of frequency bins
rng(4);
% Signal parameters
ind = 1;
while(ind)
    n = input('Select the number of source signals : ');
    if(~isnumeric(n))
        fprintf(2,'n must be an integer\n');
        ind = 1;
    else
        ind = 0;
    end
end
ind = 1;
while(ind)
    m = input('Select the number of sensors : ');
    if(~isnumeric(m))
        fprintf(2,'m must be an integer\n');
        ind = 1;
    elseif(m <= n)
        fprintf(2,'m must be larger than n\n');
        ind = 1;
    else
        ind = 0;
    end
end
if(n == 1)
    f_init = 0.0;             % LFM initial frequencies of n sources
    f_end  = 0.5;             % LFM end frequencies of n sources
    ra     = 10;              % Reception angles of n sources
elseif(n == 2)
    f_init = [0.0, 0.5];      % LFM initial frequencies of n sources
    f_end  = [0.5, 0.0];      % LFM end frequencies of n sources
    ra     = [10 30];         % Reception angles of n sources
elseif(n == 3)
    f_init = [0.0, 0.5 0.25]; % LFM initial frequencies of n sources
    f_end  = [0.5, 0.0 0.1];  % LFM end frequencies of n sources
    ra     = [10 30 50];      % Reception angles of n sources
elseif(n > 3)
    f_init = [0.0, 0.5 0.25 (fs/2).*rand(1,n-3)];  % LFM initial frequencies of n sources
    f_end  = [0.5, 0.0 0.1 (fs/2).*rand(1,n-3)];   % LFM end frequencies of n sources
    ra     = [10 30 50 randi(90,1,n-3)];           % Reception angles of n sources
end
% Reception parameters
SNR     = (-15:1:10);         % Signal-to-noise Ratio
lamda   = 150;                % Wavelength
d_space = lamda/2;            % Element spacing in the array of m antennas
rng(1);                       % Seed initialization of random noise
% Detection Parameters
perc    = 0.4;                % Percentage of the STFD maximum power to select high-energy (t,f) points
iter_N  = 2000;
theta_N = 3e2;
theta   = linspace(0, 20, theta_N);

%% Computed variables
t  = 0:1/fs:(N-1)/fs;       % Time array
f  = 0:fs/(2*M-1):fs/2;     % Frequency array

%% Signal Generation
S = signal_model(n, f_end, f_init, N, fs); % n LFM signal generation

%% Channel Mixing Model
AIC_est    = zeros(1, length(SNR));
MDL_est    = zeros(1, length(SNR));
A = inst_model(n, m, ra, lamda, d_space); % Instantaneous Mixing Model

%% Main
for ii = 1:length(SNR)
    disp(['SNR = ' num2str(SNR(ii)), ' dB'])
    for jj = 1:iter_N
        X = awgn(A*S, SNR(ii)); % Additive White Gaussian Noise for provided SNR
        % AIC Test
        AIC_est(ii) =  AIC_est(ii) + aictest(X);
        % MDL Test
        MDL_est(ii) =  MDL_est(ii) + mdltest(X);
    end
    
end
AIC_est  = AIC_est/iter_N;
MDL_est  = MDL_est/iter_N;
disp('Finished');

%% Plotting
figure;
plot(SNR,AIC_est,'bo--'); hold on;
plot(SNR,MDL_est,'rv-.'); grid on;
axis([SNR(1) SNR(end) 0.7 n+0.3])
legend('Akaike Information Criterion (AIC)','Minimum Description Length (MDL)','Location','SouthEast');
xlabel('SNR (dB)'); ylabel('Number of estimated sources')
set(gca,'fontweight','bold','fontsize',12);

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(['DOA_MTFD_Subspace_n' num2str(n)],'-depsc');
else
    return
end