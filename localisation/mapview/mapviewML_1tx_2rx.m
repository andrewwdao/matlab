%% This script creates a figure with two subplots:
% Left: 3D visualization of L(x,y) using Map3D.
% Right: Map view (Tx, Rx positions with connecting rays) using Map2D.
% The script also prints the peak coordinate and the SNR values.

clear; clc; close all;

%% User Inputs and Configurations
ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 0.1;                   % Angle resolution (degrees)
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements
OPT_GRID_DENSITY = 5; % Define a coarse grid for initial guesses
% Transmitter, receiver positions angles
area_size = 100;
pos_tx = [50, 50];                  % Tx at center
% --- Fixed rx and aoa - Scenario 1
% pos_rx = [15.98, 54.27; 67.64, 69.28]; % Two Rx positions
% aoa_act = [-6.8; 45.6];            % True AoA from Rx to Tx
% aoa_act = [-1.1; -5.4];
% --- Fixed rx and aoa - optimal scenario
% pos_rx = [21, 51; 21, 60]; % 1
% pos_rx = [21, 51; 30, 70]; % 2
% pos_rx = [21, 51; 40, 70]; % 3
% pos_rx = [21, 51; 50, 70]; % 4
% pos_rx = [21, 51; 60, 70]; % 5
% pos_rx = [21, 51; 70, 60]; % 6
% pos_rx = [21, 51; 70, 50]; % 7
% pos_rx = [21, 51; 70, 40]; % 8
% pos_rx = [21, 51; 60, 30]; % 9
% pos_rx = [21, 51; 50, 30]; % 10
% pos_rx = [21, 51; 40, 30]; % 11
% pos_rx = [21, 51; 30, 30]; % 12
% % pos_rx = [21, 51; 0, 30]; % 13
pos_rx = [21, 51; 20, 40; 70,50; 60,80]; % 13
aoa_act = [0; 0; 0; 0];            % True AoA from Rx to Tx
RX_NUM = length(pos_rx);                         % Number of receivers

% --- Randomised rx and aoa
% RX_NUM = 5;
% pos_rx = area_size*rand(RX_NUM, 2); % Random Receiver position (x, y) in meters
% aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1);

% Physical constants and wavelength
SNR_dB = 10 * ones(RX_NUM, 1);      % SNR in dB
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength
progressbar('reset', RX_NUM+1); % Reset progress bar
progressbar('minimalupdateinterval', 0); % Set a smaller interval at the beginning
SHOW_LIMITS = false; % Show the detecting limits of the RXs (with known limitation)
SHOW_EXTRA = true; % Show extra information such as the AoA and the intersection point

%% Compute absolute angles from each Rx to Tx and corresponding rotations
angle_rx_tx_abs = zeros(RX_NUM, 1);
for i = 1:RX_NUM
    angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
end
rot_abs = angle_rx_tx_abs - aoa_act;

%% Signal and channel configurations
avg_amp_gain = 1;
P_t = ones(RX_NUM, 1);
sub_carrier = (1:RX_NUM)' * 1000;
Fs = 2 * max(sub_carrier);
T = TIME_INST_NUM / Fs;
t = 0:1/Fs:(T-1/Fs);
element_spacing = 0.5 * lambda;
sweeping_angle = -90:RESOLUTION:90;

%% For each Rx, generate received signal
w = cell(RX_NUM, 1);
% Initialize channel model
channel = ChannelModels();
for rx_idx = 1:RX_NUM
    % Generate base signal (using last element index as in the original code)
    s_t = sqrt(P_t(RX_NUM)) .* exp(1j * 2 * pi * sub_carrier(RX_NUM) * t);
    avg_E = FIXED_TRANS_ENERGY * 1 + ~FIXED_TRANS_ENERGY * (avg_amp_gain^2 * P_t(RX_NUM) * T * Fs);
    nPower = avg_E / db2pow(SNR_dB(rx_idx));
    y_los = channel.LoS(s_t, avg_amp_gain);
    y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
    y_awgn = channel.AWGN(y_ula, nPower);
    w{rx_idx} = y_awgn;
    progressbar('step'); % Update progress bar
end
%% Find the Maximum Likelihood (ML) estimate of the transmitter position
% nPower_model = 1; % Noise power level for the model
l4c = Likelihood4Coordinates();
optimiser = gridOptimiser();
objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
[optCoord, L_peak] = optimiser.fmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
progressbar('end');  % This will display the total runtime
progressbar('reset', 1); % Reset progress bar
%% additional function to plot the 3D map
[X, Y, L] = l4c.calLikelihood4Area(area_size, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
progressbar('end');  % This will display the total runtime

%% === Plotting and display the result
fprintf('Peak found at (%.2f, %.2f) with L = %.2f\n', optCoord(1), optCoord(2), L_peak);
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
map2d.plot(pos_tx, pos_rx, rot_abs, area_size, aoa_act, ABS_ANGLE_LIM, [SHOW_LIMITS, SHOW_EXTRA]);
