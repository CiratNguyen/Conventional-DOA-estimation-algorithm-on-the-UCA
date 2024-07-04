close all
clear
clc
%% (1) Initializing signal model, UCA configuration and simulation parameters
Fs = 48000; Ts = 1/Fs; % Sampling frequency (Hz), Sampling period (s)
ts = Ts:Ts:0.1; % Sampling time (s)
f1 = 1325; s_A = 1*sin(2*pi*f1*ts).'; % Signal A's frequency (Hz), Signal A's amplitude (mV)
f2 = 1350; s_B = 1*sin(2*pi*f2*ts).'; % Signal B's frequency (Hz), Signal B's amplitude (mV)
f3 = 1375; s_C = 1*sin(2*pi*f3*ts).'; % Signal C's frequency (Hz), Signal C's amplitude (mV)
Sn = length(ts); % Snapshots
SNR = 30; % Signal to noise ratio (SNR)
angles_A = -120*pi/180; % Signal A's DOA [azimuth]
angles_B = 0*pi/180; % Signal B's DOA [azimuth]
angles_C = 70*pi/180; % Signal C's DOA [azimuth]
c = 343; % Sound propagation's speed in the air 
lambda_max = c/min([f1 f2 f3]); % Lambda max
Ne = 16; % Number of elements
D = 3; % Number of signals
R = 0.045; % Array's radius
k = 2*pi/lambda_max; %  Array's angular coefficient
k1 = 2*pi/(c/f1); %  Signal A's angular coefficient
k2 = 2*pi/(c/f2); %  Signal B's angular coefficient
k3 = 2*pi/(c/f3); %  Signal C's angular coefficient
phi=-180:0.05:180;% Azimuth plane
for l = 1:Ne
    a_A(l)=exp(1j*k1*R*cos(angles_A-2*pi*((l-1)/Ne))); %Signal A's steering vector at array
    a_B(l)=exp(1j*k2*R*cos(angles_B-2*pi*((l-1)/Ne))); %Signal B's steering vector at array
    a_C(l)=exp(1j*k3*R*cos(angles_C-2*pi*((l-1)/Ne))); %Signal C's steering vector at array
end
for l = 1:Ne
    x_A(l,:)=a_A(l)*s_A; % Signal A at array
    x_B(l,:)=a_B(l)*s_B; % Signal B at array
    x_C(l,:)=a_C(l)*s_C; % Signal C at array
end
x = x_A+x_B+x_C; %x0=∑x0i (i=1,2,...,Ne); Tổng tín hiệu tại mảng
x = awgn(x,SNR,'measured'); % Thêm nhiễu Gaussian trắng vào tập tín hiệu
%% (2) Decompositing the covariance matrix
Rx = x*x'/Sn;  % Calculating the covariance matrix
[eigvec,eigval] = eig(Rx); % Tính giá trị riêng và vector riêng
En = eigvec(:,1:Ne-D); % Decompositing the eigenvalues and the eigenvectors
%% (3) Using MUSIC estimation algorithm for estimating DOA of signal(s)
for pp = 1:length(phi)
    for l = 1:Ne
        a0(l,1) = exp(1j*k*R*cos(phi(pp)*pi/180-2*pi*((l-1)/Ne)));
    end
    Pmusic(1,pp) = abs(1/(a0'*En*(En')*a0)); 
end
%% (4) Estimating DOA of signal(s) by MUSIC
figure(1);
plot(phi,10*log10(Pmusic/max(Pmusic))); title('MUSIC algorithm');
xlabel('Azimuth plane'); ylabel('Amplitude (dB)'); grid on;