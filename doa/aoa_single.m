clear; clc; close all;
%% ====================== User inputs
ITERATION = 100;
TRUE_ANGLE = 30;
SNR_db = 10; %dB

%% ====================== Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
t = (0:1e-6:1e-3);  % Time vector for the signal
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
tx_pos = [10+10*cosd(TRUE_ANGLE), 50+10*sind(TRUE_ANGLE);]; % Transmitter position (x, y) in meters
rx_pos = [10, 50;]; % Receiver position (x, y) in meters
element_num = 4;   % Number of elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:1:90; % Angle range for finding the AoA
P_t = 1; % W - Transmit signal power
nPower = P_t/db2pow(SNR_db); % W - White noise power
%% ====================== Generate original signal received at Rx
% re = randn; im = randn;
% s_t = sqrt(P_t) * (re+1j*im) / sqrt(re^2+im^2); 
s_t = sqrt(P_t); % Generate transmitted signal
channel = ChannelModel(tx_pos, rx_pos, lambda, element_num, lambda/2);
array = ArrayModel(tx_pos, rx_pos, lambda, element_num, lambda/2);
y_t_original = channel.Direct(s_t);  % Received signal at the receiver
y_t_original = array.applyULA(y_t_original);  % Apply ULA characteristics to the received signal
%% ====================== Apply noise
noise = sqrt(nPower/2) * (randn(size(y_t_original)) + 1j * randn(size(y_t_original))); % White noise
y_t = y_t_original + noise;

%% ---------------------- DoA Estimation Algorithm ----------------------
estimator = DoAEstimator(y_t, size(tx_pos,1), lambda, element_num, element_spacing, sweeping_angle);
[est_aoa_sync, spectrum_dB_sync] = estimator.ML_sync();
[est_aoa_bf, spectrum_dB_bf] = estimator.BF();
[est_aoa_mvdr, spectrum_dB_mvdr] = estimator.MVDR();
[est_aoa_music, spectrum_dB_music] = estimator.MUSIC();

%% --- Plotting
figure; hold on; grid on;
plot(spectrum_dB_sync);
plot(spectrum_dB_bf);
vis_sync = DoAVisualisation(...
    "ML Sync", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB_sync, ...
        est_aoa_sync);
vis_sync.plot();
vis_bf = DoAVisualisation(...
    "Beamforming", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB_bf, ...
        est_aoa_bf);
vis_bf.plot();
vis_mvdr = DoAVisualisation(...
    "MVDR", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB_mvdr, ...
        est_aoa_mvdr);
vis_mvdr.plot();
vis_music = DoAVisualisation(...
    "MUSIC", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB_music, ...
        est_aoa_music);
vis_music.plot();
