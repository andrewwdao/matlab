%Initialisation and visualise Tx/Rx on a map
clear; clc; close all;
rs=rng(2007); % initialize the random number generator to a specific seed value
% Constants
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength

% Area size and positions
area_size = 100;   % 100x100 meter area
tx_pos = [60,80;]; % Transmitter position (x, y) in meters
rx_pos = [10,50;]; % Receiver position (x, y) in meters
element_num = 4;   % Number of elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:1:90; % Angle range for sweeping to find the AoA

t = (0:1e-6:1e-3);  % Time vector for the signal
P_t = 1;  % Transmit signal power
% Signal and noise parameters
Nsamp = 1;
nPower_db = 10; % White noise power (dB)

% ====================== Generate signal received at Rx =======================
% --- Generate transmitted signal
% Each column of X represents the received signal at the corresponding element.
% Transmitted signal (simple sinusoidal signal)
% s_t = sqrt(P_t) * exp(1j * 2 * pi * fc * t);  % Complex sinusoid
s_t = sqrt(P_t) * (1+1j); % --- Considering one time instance only
% --- Transmitted through a Friis free-space model
% Calculate the Euclidean distance
channel = ChannelModelULA(tx_pos, rx_pos, lambda, element_num, lambda/2, 0);
y_t = channel.FriisModel(s_t);  % Received signal at the receiver


%% ---------------------- DoA Estimation Algorithm ----------------------
estimator = DoAEstimator(y_t, size(tx_pos,1), lambda, element_num, element_spacing, sweeping_angle, db2pow(nPower_db));
[est_aoa, spectrum_dB] = estimator.MUSIC();


% % --- Plotting
vis = DoAVisualisation(...
    "MUSIC", tx_pos, rx_pos, area_size, ...
    sweeping_angle, spectrum_dB, ...
        est_aoa);
vis.plot();
