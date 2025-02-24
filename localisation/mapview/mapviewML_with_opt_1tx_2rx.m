%% mapview_with_visualisation.m
% This script creates a figure with two subplots:
% Left: 3D visualization of L(x,y) using Map3D.
% Right: Map view (Tx, Rx positions with connecting rays) similar to mapview_1tx_2rx.m.

clear; clc; close all;

%% User Inputs and Configurations
RX_NUM = 2;                         % Number of receivers
SNR_dB = 10 * ones(RX_NUM, 1);       % SNR in dB
ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 0.1;                   % Angle resolution (degrees)
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements

% Physical constants and wavelength
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength

% Transmitter, receiver positions angles
area_size = 100;
pos_tx = [50, 50];                  % Tx at center
% --- Randomised rx and aoa
% pos_rx = area_size*rand(RX_NUM, 2); % Random Receiver position (x, y) in meters
% aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1);
% --- Fixed rx and aoa - Scenario 1
% pos_rx = [15.98, 54.27; 67.64, 69.28]; % Two Rx positions
% aoa_act = [-6.8; 45.6];            % True AoA from Rx to Tx
% aoa_act = [-1.1; -5.4];
% --- Fixed rx and aoa - optimal scenario
pos_rx = [21, 51; 51, 81]; % Two Rx positions
aoa_act = [0; 0];            % True AoA from Rx to Tx


% Compute absolute angles from each Rx to Tx and corresponding rotations
angle_rx_tx_abs = zeros(RX_NUM, 1);
for i = 1:RX_NUM
    angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
end
rot_abs = angle_rx_tx_abs - aoa_act;

% Signal and channel configurations
avg_amp_gain = 1;
P_t = ones(RX_NUM, 1);
sub_carrier = (1:RX_NUM)' * 1000;
Fs = 2 * max(sub_carrier);
T = TIME_INST_NUM / Fs;
t = 0:1/Fs:(T-1/Fs);
element_spacing = 0.5 * lambda;
sweeping_angle = -90:RESOLUTION:90;

% For each Rx, generate received signal (only channel initialization is shown)
w = cell(RX_NUM, 1);
% Initialize channel model
channel = ChannelModels();
for rx_idx = 1:RX_NUM
    % Generate base signal (using last element index as in the original code)
    s_t = sqrt(P_t(RX_NUM)) .* exp(1j * 2 * pi * sub_carrier(RX_NUM) * t);
    avg_E = FIXED_TRANS_ENERGY * 1 + ~FIXED_TRANS_ENERGY * (avg_amp_gain^2 * P_t(RX_NUM) * T * Fs);
    % nPower = avg_E / db2pow(SNR_dB(rx_idx));
    nPower = 0;
    y_los = channel.LoS(s_t, avg_amp_gain);
    y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
    y_awgn = channel.AWGN(y_ula, nPower);
    w{rx_idx} = y_awgn;
end

nPower_model = 1; % Noise power level for the model
l4c = Likelihood4Coordinates();
[X, Y, L] = l4c.CalculateLikelihood4Area(area_size, pos_rx, rot_abs, w, ELEMENT_NUM, nPower_model);
% Define the objective function to maximize
objective_to_maximize = @(coor) -abs(l4c.fminconCalculateLikelihood(coor, pos_rx, rot_abs, w, ELEMENT_NUM, nPower_model));

% Set bounds
lb = [0, 0];
ub = [area_size, area_size];

% Use MLoptimiser class to find the maximum likelihood estimate
grid_points = 5; % Define a 5x5 coarse grid for initial guesses
optimizer = gridOptimiser();
[optCoord, L_peak] = optimizer.fmincon(objective_to_maximize, {}, lb, ub, grid_points);

% Print the result
fprintf('Peak found at (%.2f, %.2f) with L = %.2f\n', optCoord(1), optCoord(2), L_peak);

%% === Plotting
fprintf('SNR (dB):\n')
fprintf('%.0f  ', SNR_dB); fprintf('\n');
fprintf('RX Positions:\n');
for idx = 1:size(pos_rx, 1)
    fprintf('  Rx %d: x = %.2f, y = %.2f\n', idx, pos_rx(idx, 1), pos_rx(idx, 2));
end
fprintf('Absolute RX Rotations (degrees):\n  ');
fprintf('%.2f  ', rot_abs); fprintf('\n');
fprintf('Absolute RX-TX Angles (degrees):\n  ');
fprintf('%.2f  ', angle_rx_tx_abs); fprintf('\n');
fprintf('True AoA (degrees):\n  ');
fprintf('%.2f  ', aoa_act); fprintf('\n');

figure('Name', '3D ML Visualization and Map View', 'WindowState', 'maximized');
% Left Subplot: 3D Visualization of ML Function
subplot(1,2,1);
map3d = Map3D(X, Y, L);
map3d.plot(gca);
% Right Subplot: Map View
subplot(1,2,2); hold on;
map2d = Map2D();
map2d.plot(pos_tx, pos_rx, rot_abs, area_size, aoa_act, ABS_ANGLE_LIM, false);

