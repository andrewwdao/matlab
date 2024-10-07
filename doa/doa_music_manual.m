%Initialisation and visualise Tx/Rx on a map
clear; clc; close all;
rs=rng(2007); % initialize the random number generator to a specific seed value
% Constants
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength

% Area size and positions
area_size = 100;   % 100x100 meter area
tx_pos = [60,80; 40,40;]; % Transmitter position (x, y) in meters
rx_pos = [10,50;]; % Receiver position (x, y) in meters
element_num = 4;   % Number of elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:1:90; % Angle range for sweeping to find the AoA

t = (0:1e-6:1e-3);  % Time vector for the signal
P_t = 1;  % Transmit signal power
nPower_db = -90; % White noise power (dB)

% ====================== Generate signal received at Rx =======================
% --- Generate transmitted signal
% Each column of X represents the received signal at the corresponding element.
% s_t = sqrt(P_t) * exp(1j * 2 * pi * fc * t);  % Complex sinusoid
s_t = sqrt(P_t) * (1+1j); % --- Considering one time instance only
% --- Transmitted through a Friis free-space model with ULA characteristics at the receiver
channel = ChannelModel(tx_pos, rx_pos, lambda, element_num, lambda/2, 0);
y_t = channel.FriisModel(s_t);  % Received signal at the receiver
array = ArrayModel(tx_pos, rx_pos, lambda, element_num, lambda/2, 0);
y_t = array.applyULA(y_t);  % Apply ULA characteristics to the received signal
noise = sqrt(db2pow(nPower_db)/2) * (randn(size(y_t)) + 1j * randn(size(y_t))); % White noise
y_t = y_t + noise;

%% ---------------------- DoA Estimation Algorithm ----------------------
estimator = DoAEstimator(y_t, size(tx_pos,1), lambda, element_num, element_spacing, sweeping_angle, db2pow(nPower_db));
[est_aoa_music, spectrum_dB_music] = estimator.MUSIC();
[est_aoa_mvdr, spectrum_dB_mvdr] = estimator.MVDR();
[est_aoa_beam, spectrum_dB_beam] = estimator.Beamforming();

%% --- Plotting
vis_music = DoAVisualisation(...
    "MUSIC", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB_music, ...
        est_aoa_music);
vis_music.plot();
vis_mvdr = DoAVisualisation(...
    "MVDR", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB_mvdr, ...
        est_aoa_mvdr);
vis_mvdr.plot();
vis_beam = DoAVisualisation(...
    "Beamforming", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB_beam, ...
        est_aoa_beam);
vis_beam.plot();
