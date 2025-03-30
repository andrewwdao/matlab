%% Performs Maximum Likelihood (ML) localization of a single randomly
% positioned transmitter using a fixed array of receivers with Uniform Linear Arrays (ULAs).
%

% Overview:
% This script simulates localization of a randomly positioned transmitter using multiple
% fixed receivers equipped with Uniform Linear Arrays (ULAs). It implements Maximum Likelihood (ML)
% estimation to determine the transmitter's position based on Angle of Arrival (AoA) measurements.
%
% Key Components:
% 1. Environment Setup: Creates a 2D environment with fixed receivers and a random transmitter
%    position with safety distance constraints
% 2. Signal Generation: Models ULA signal reception with configurable SNR, carrier frequency,
%    and element spacing
% 3. ML Estimation: Computes likelihood function across the area and finds optimal transmitter
%    position using grid-based optimization
% 4. Visualization: Displays 3D likelihood surface and 2D map showing receiver positions,
%    transmitter location, and AoA information
%
% Dependencies:
%   - Map2D class: Handles 2D visualization and position generation
%   - Map3D class: Provides 3D visualization of likelihood function
%   - ChannelModels class: Implements signal propagation models
%   - Likelihood4Coordinates class: Calculates likelihood for position estimation
%   - Optimisers class: Provides optimization algorithms for ML estimation
%   - progressbar function: Displays progress during computation
%
% Outputs:
%   - Console output: True TX position, estimated position, positioning error,
%     RX positions, absolute RX rotations, and true AoA values
%   - Figure with two subplots:
%     * Left: 3D visualization of ML likelihood function L(x,y)
%     * Right: 2D map view showing receiver positions, true transmitter position,
%       AoA measurements, and detection limits (optional)
%
% Author: Minh An Dao (Andreww)

clear; clc; close all;

%% User Inputs and Configurations
ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 0.1;                   % Angle resolution (degrees)
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements
OPT_GRID_DENSITY = 5;               % Define a coarse grid for initial guesses
RX_NUM = 10;                        % Number of receivers
SAFETY_DISTANCE = 2;                % Minimum distance between TX and RX (meters)

% Area definition
area_size = 100;
x_min = 0;
x_max = area_size;
y_min = 0;
y_max = area_size;

% Physical constants and wavelength
SNR_dB = 50 * ones(RX_NUM, 1);      % SNR in dB
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength

progressbar('reset', RX_NUM+1);     % Reset progress bar
progressbar('minimalupdateinterval', 0); % Set a smaller interval at the beginning

SHOW_LIMITS = false;                % Show the detecting limits of the RXs
SHOW_EXTRA = true;                  % Show extra information such as AoA and intersection

%% Initialise Map2D with fixed receivers
% Fixed receiver positions between [10,10] and [90,90]
map2d = Map2D([10,10], [90, 90], RX_NUM);

%% Generate random transmitter position
% Random position within the area bounds with safety margin
margin = SAFETY_DISTANCE;
pos_tx = [rand(1)*(x_max-x_min-2*margin)+x_min+margin, ...
          rand(1)*(y_max-y_min-2*margin)+y_min+margin];

fprintf('Random transmitter position: (%.2f, %.2f)\n', pos_tx(1), pos_tx(2));

%% Compute absolute angles from each Rx to Tx and corresponding rotations
[pos_rx, aoa_act, rot_abs] = map2d.genPos(area_size, pos_tx, RX_NUM, false, SAFETY_DISTANCE, ABS_ANGLE_LIM, RESOLUTION);

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
    % Generate base signal
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
l4c = Likelihood4Coordinates();
optimiser = Optimisers();
objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
[optCoord, L_peak] = optimiser.gridFmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
progressbar('end');  % This will display the total runtime
progressbar('reset', 1); % Reset progress bar

%% Additional function to plot the 3D map
[X, Y, L] = l4c.calLikelihood4Area(area_size, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
progressbar('end');  % This will display the total runtime

%% === Plotting and display the result
fprintf('True TX position: (%.2f, %.2f)\n', pos_tx(1), pos_tx(2));
fprintf('Estimated TX position: (%.2f, %.2f) with L = %.2f\n', optCoord(1), optCoord(2), L_peak);
fprintf('Positioning error: %.2f meters\n', norm(pos_tx - optCoord));

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