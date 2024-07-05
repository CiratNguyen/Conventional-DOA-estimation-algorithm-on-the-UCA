close all
clear
clc

%% (1) Initializing signal model, UCA configuration and simulation parameters
Fs = 48000; Ts = 1/Fs; % Sampling frequency (Hz), Sampling period (s)
ts = Ts:Ts:0.1; % Sampling time (s)
f = 1300; s = 1*sin(2*pi*f*ts).'; % Signal's frequency (Hz), Signal's amplitude (mV)
Sn = length(ts); % Snapshots
DOA = 0*pi/180; % Signal's DOA
c = 343; % Sound propagation's speed in the air
lambda_max = c/f; % Lambda max
Ne = 6; % Number of elements
D = 1; % Number of signal(s)
R = 0.5*Ne*lambda_max/(2*pi); % Array's radius
ka = 2*pi/lambda_max; %  Array's angular coefficient
ks = 2*pi/(c/f); %  Signal's angular coefficient
phi=-180:0.05:180; % Azimuth plane
for l = 1:Ne
    a(l)=exp(1j*ks*R*cos(DOA-2*pi*((l-1)/Ne))); %Signal's steering vector at array
end
for l = 1:Ne
    x(l,:)=a(l)*s; %Signal's at array
end
SNRs = -10:1:30; % Signal to noise ratio (SNR)
for SNR_idx = 1:length(SNRs)
    SNR = SNRs(SNR_idx);
    for ite = 1:1
        x1 = x;
        x2 = awgn(x1,SNR,'measured'); % Adding White Gaussian noise to signal
        %% (2) Decompositing the covariance matrix
        Rx = x2*x2'/Sn; % Calculating the covariance matrix
        [eigvec,eigval]=eig(Rx); % Decompositing the eigenvalues and the eigenvectors
        En = eigvec(:,1:Ne-D); % Constructing the noise subspace
        invRx = inv(Rx); %InverseRx for MVDR
        %% (3) Using DOA estimation algorithms for estimating DOA of signal(s)
        for pp = 1:length(phi)
            for l = 1:Ne
                a1(l,1) = exp(1j*ka*R*cos(phi(pp)*pi/180-2*pi*((l-1)/Ne)));
            end
            Pmusic(1,pp) = abs(1/(a1'*En*(En')*a1)); %MUSIC algorithm spatial spectrum
            Pcb(1,pp) = abs(a1'*Rx*a1); %CB algorithm spatial spectrum
            Pmvdr(1,pp) = abs(1/(a1'*invRx*a1)); %MVDR algorithm spatial spectrum
        end
        %% (4) Estimating DOA of signal(s) by MUSIC, CB and MVDR
        [u1,v1] = max(Pmusic);
        DOA_music = phi(v1);
        PAPR_music(ite) = u1/mean(Pmusic);
        RMSE_music(ite) = DOA - DOA_music;
        [u2,v2] = max(Pcb);
        DOA_cb = phi(v2);
        PAPR_cb(ite) = u2/mean(Pcb);
        RMSE_cb(ite) = DOA - DOA_cb;
        [u3,v3] = max(Pmvdr);
        DOA_mvdr = phi(v3);
        PAPR_mvdr(ite) = u3/mean(Pmvdr);
        RMSE_mvdr(ite) = DOA - DOA_mvdr;
    end
    mean_PAPR(1,SNR_idx) = 10*log10(mean(PAPR_music));
    mean_PAPR(2,SNR_idx) = 10*log10(mean(PAPR_cb));
    mean_PAPR(3,SNR_idx) = 10*log10(mean(PAPR_mvdr));
    DOA_RMSE(1,SNR_idx) = sqrt(mean(RMSE_music.^2));
    DOA_RMSE(2,SNR_idx) = sqrt(mean(RMSE_cb.^2));
    DOA_RMSE(3,SNR_idx) = sqrt(mean(RMSE_mvdr.^2));
end
figure(1);
plot(SNRs,mean_PAPR(1,:)); hold on; plot(SNRs,mean_PAPR(2,:)); hold on; plot(SNRs,mean_PAPR(3,:)); grid on;
legend('MUSIC','CB','MVDR');
figure(2);
plot(SNRs,DOA_RMSE(1,:)); hold on; plot(SNRs,DOA_RMSE(2,:));hold on; plot(SNRs,DOA_RMSE(3,:));
legend('CB','MUSIC','MVDR'); grid on;
