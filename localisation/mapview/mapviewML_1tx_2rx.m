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
RX_NUM = 10;                        % Number of receivers
RANDOMISE_RX = true;                % Randomise RX positions and AoA
SAFETY_DISTANCE = 2;                % Minimum distance between TX and RX (meters)
% Transmitter, receiver positions angles
area_size = 100;
pos_tx = [50, 50];                  % Tx at center of area
% Physical constants and wavelength
SNR_dB = 50 * ones(RX_NUM, 1);      % SNR in dB
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength
progressbar('reset', RX_NUM+1); % Reset progress bar
progressbar('minimalupdateinterval', 0); % Set a smaller interval at the beginning
SHOW_LIMITS = false; % Show the detecting limits of the RXs (with known limitation)
SHOW_EXTRA = true; % Show extra information such as the AoA and the intersection point
%% Initialise classes
map2d = Map2D();
%% Compute absolute angles from each Rx to Tx and corresponding rotations
[pos_rx, aoa_act, rot_abs] = map2d.genPos(area_size, pos_tx, RX_NUM, RANDOMISE_RX, SAFETY_DISTANCE, ABS_ANGLE_LIM, RESOLUTION);
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
optimiser = Optimisers();
objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
[optCoord, L_peak] = optimiser.gridFmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
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
fprintf('True AoA (degrees):\n  ');
fprintf('%.2f  ', aoa_act); fprintf('\n');

figure('Name', '3D ML Visualization and Map View', 'WindowState', 'maximized');
% Left Subplot: 3D Visualization of ML Function
subplot(1,2,1);
map3d = Map3D(X, Y, L);
map3d.plot(gca);
% Right Subplot: Map View
subplot(1,2,2); hold on;
map2d.plot(pos_tx, pos_rx, rot_abs, area_size, aoa_act, ABS_ANGLE_LIM, [SHOW_LIMITS, SHOW_EXTRA]);
