clear; clc; close all;
%% ====================== User inputs
ITERATION = 1;
TRUE_ANGLE = 30;
SNR_dB = 10; %dB
avg_amp_gain = 1; % Average gain of the channel

%% ====================== Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
t = (0:1e-6:1e-3);  % Time vector for the signal
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
tx_pos = [10+10*cosd(TRUE_ANGLE), 50+10*sind(TRUE_ANGLE);]; % Transmitter position (x, y) in meters
P_t = 1;  % W - Transmit signal power1000Hz
rx_pos = [10, 50;]; % Receiver position (x, y) in meters
element_num = 4;   % Number of elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:1:90; % Angle range for finding the AoA
%% ====================== Generate original signal received at Rx
s_t = sqrt(P_t(1)) .* exp(1j * 2 * pi);
channel = ChannelModel(tx_pos, rx_pos, lambda, element_num, lambda/2);
y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
y_ula = channel.applyULA(y_los);  % Apply ULA characteristics to the received signal
% Calculate the energy of the whole signal transmitted to correct the noise power
avg_E =  avg_amp_gain^2 * P_t; % average received signal energy, \sum_{t=1}^{T} ||s(t)||^2
nPower = avg_E/db2pow(SNR_dB); % W - White noise power - noise variance
y_awgn = channel.AWGN(y_ula, nPower);

%% ---------------------- DoA Estimation Algorithm ----------------------
estimator = DoAEstimator(y_awgn, size(tx_pos,1), lambda, element_num, element_spacing, sweeping_angle);
[est_aoa_sync, spectrum_dB_sync] = estimator.ML_sync(s_t);
[est_aoa_bf, spectrum_dB_bf] = estimator.BF();
[est_aoa_mvdr, spectrum_dB_mvdr] = estimator.MVDR();
[est_aoa_music, spectrum_dB_music] = estimator.MUSIC();

%% --- Plotting
figure; hold on; grid on;
plot(spectrum_dB_sync, 'LineWidth', 1, 'DisplayName', "Sync");
plot(spectrum_dB_bf, 'LineWidth', 1, 'DisplayName', "BF");
title("Spectrum"); legend("AutoUpdate","on");
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
