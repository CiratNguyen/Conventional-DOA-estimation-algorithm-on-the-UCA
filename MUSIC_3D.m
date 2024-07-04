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
angles_A = [-120 10]*pi/180; % Signal A's DOA [azimuth elevation]
angles_B = [0 40]*pi/180; % Signal B's DOA [azimuth elevation]
angles_C = [70 50]*pi/180; % Signal C's DOA [azimuth elevation]
c = 343; % Sound propagation's speed in the air
lambda_max = c/min([f1 f2 f3]); % Lambda max
Ne = 6; % Number of elements
D = 3;  % Number of signals
R = 0.5*Ne*lambda_max/(2*pi); % Array's radius
k = 2*pi/lambda_max; %  Array's angular coefficient
k1 = 2*pi/(c/f1); %  Signal A's angular coefficient
k2 = 2*pi/(c/f2); %  Signal B's angular coefficient
k3 = 2*pi/(c/f3); %  Signal C's angular coefficient
phi = -180:0.05:180; % Azimuth plane
theta = 0:0.05:90; % Elevation plane
for l = 1:Ne
    a_A(l) = exp(1j*k1*R*sin(angles_A(2))*cos(angles_A(1)-2*pi*((l-1)/Ne))); %Signal A's steering vector at array
    a_B(l) = exp(1j*k2*R*sin(angles_B(2))*cos(angles_B(1)-2*pi*((l-1)/Ne))); %Signal B's steering vector at array
    a_C(l) = exp(1j*k3*R*sin(angles_C(2))*cos(angles_C(1)-2*pi*((l-1)/Ne))); %Signal C's steering vector at array
end
for l = 1:Ne
    x_A(l,:) = a_A(l)*s_A; % x0_A=a0_A*s0_A tín hiệu A tại mảng
    x_B(l,:) = a_B(l)*s_B; % x0_B=a0_B*s0_B tín hiệu B tại mảng
    x_C(l,:) = a_C(l)*s_C; % x0_C=a0_C*s0_C tín hiệu C tại mảng
end
x = x_A+x_B+x_C; %x0=∑x0i (i=1,2,...,Ne); Tổng tín hiệu tại mảng
x = awgn(x,SNR,'measured'); % Thêm nhiễu Gaussian trắng vào tập tín hiệu
%% (2) Decompositing the covariance matrix
Rx = x*x'/Sn;  % Calculating the covariance matrix
[eigvec,eigval] = eig(Rx); % Tính giá trị riêng và vector riêng
En = eigvec(:,1:Ne-D); % Decompositing the eigenvalues and the eigenvectors
%% (3) Using MUSIC estimation algorithm for estimating DOA of signal(s)
for tt = 1:length(theta)
    for pp = 1:length(phi) 
        for l = 1:Ne
            a0(l,1) = exp(1j*k*R*sin(theta(tt)*pi/180)*cos(phi(pp)*pi/180-2*pi*((l-1)/Ne)));
        end
        Pmusic(tt,pp) = abs(1/(a0'*En*(En')*a0)); %cong thuc thuat toan MUSIC
    end
end
%% (4) Estimating DOA of signal(s) by MUSIC
figure(1); title('MUSIC algorithm')
surf(phi,theta,10*log10(Pmusic/max(Pmusic(:))),'EdgeColor', 'none');
colormap('jet'); colorbar; xlabel('Azimuth plane'); ylabel('Elevation plane'); grid on;
local_maxima1 = imregionalmax(Pmusic);
[row1, col1] = find(local_maxima1);
[max_values1, indices1] = maxk(Pmusic(local_maxima1), D);
max_row1 = sort(row1(indices1));
max_col1 = sort(col1(indices1));
fprintf('MUSIC algorithm\n');
for i = 1:length(max_row1)
    fprintf('DOA %d:\n [%f, %f]\n', i, phi(max_col1(i)), theta(max_row1(i)));
end




