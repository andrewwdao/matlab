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
DOA_MODE = 'sweep';                 % DoA estimation mode ('sweep' or 'opt')
DOA_RESOLUTION = 1;               % Angle resolution (degrees)
OPT_GRID_DENSITY = 5;               % Define a coarse grid for initial guesses
NUM_RX = 10;                        % Number of receivers
NUM_TX = 1;                         % Number of transmitters
RANDOMISE_TX = true;               % Randomise TX positions
SAFETY_DISTANCE = 5;                % Minimum distance between TX and RX (meters)
area_size = 100;
% Physical constants and wavelength
SNR_dB = 50 * ones(1, NUM_RX);      % SNR in dB
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength
avg_amp_gain = 1;                           % Average gain of the channel
L_d0=100;                                   % Reference Power (dB) - for gain calculation
d0=100;                                     % Reference distance (m) - for gain calculation
alpha=4;                                    % Path loss exponent - for gain calculation
P_t = 1;                                    % W - Transmit signal power (known)                      
Fs = 2 * fc;                                % Sample frequency, enough for the signal
T = TIME_INST_NUM/Fs;                       % Period of transmission
t = 0:1/Fs:(T-1/Fs);                        % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;
sweeping_angle = -90:RESOLUTION:90;

SHOW_LIMITS = false;                % Show the detecting limits of the RXs
SHOW_EXTRA = true;                  % Show extra information such as AoA and intersection

%% Initialise classes
channel = ChannelModels();
l4c = Likelihood4Coordinates();
optimiser = Optimisers();
algo = Algorithms(l4c, optimiser);
map2d = Map2D([10,10], [90, 90], NUM_RX);
map3d = Map3D(map2d);
%% Generate random transmitter position
pos_tx = map2d.genTXPos(area_size, NUM_TX, RANDOMISE_TX);
fprintf('Random transmitter position: (%.2f, %.2f)\n', pos_tx(1), pos_tx(2));
%% Compute absolute angles from each Rx to Tx and corresponding rotations
[pos_rx, aoa_act, rot_abs] = map2d.genRXPos(area_size, pos_tx, NUM_RX, false, SAFETY_DISTANCE, ABS_ANGLE_LIM, RESOLUTION);
% Generate nuisance transmitted signal with random phase
[s_t, e_avg] = channel.generateNuisanceSignal(fc, P_t, T, t, TIME_INST_NUM, FIXED_TRANS_ENERGY);
[nPower, y_centralised] = channel.generateReceivedSignal(s_t, pos_tx, pos_rx, aoa_act, e_avg, SNR_dB, L_d0, d0, alpha, ELEMENT_NUM, element_spacing, lambda);
%% Initialise arrays
ula = ULA(lambda, ELEMENT_NUM, element_spacing);    % Create Uniform Linear Array object
estimator = DoAEstimator(ula, sweeping_angle, 0, DOA_MODE, OPT_GRID_DENSITY);
doa_estimator = @(sig) estimator.BF(sig);

%% Alternative display for the RXs - one receiver per line
fprintf('Receiver Information (fixed):\n');
fprintf('%-6s %-15s %-15s %-15s\n', 'RX', 'Position', 'Rotation (°)', 'AoA (°)');
fprintf('----------------------------------------------------------------\n');

for idx = 1:size(pos_rx, 1)
    fprintf('%-6d (%-5.2f, %-5.2f) %-15.2f %-15.2f\n', ...
        idx, pos_rx(idx,1), pos_rx(idx,2), rot_abs(idx), aoa_act(idx));
end
fprintf('\n');

%% Find the Maximum Likelihood (ML) estimate of the transmitter position
progressbar('reset', 1);     % Reset progress bar
[optCoord, L_peak, error] = algo.MLOpt4mCentroid(...
    pos_rx, rot_abs, y_centralised(1,:,:), ...
    ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
    doa_estimator, pos_tx...
);
progressbar('end');  % This will display the total runtime

% Add TX and estimation info
annotStrings = {};
annotStrings{end+1} = sprintf('True TX: (%.2f, %.2f)', pos_tx(1), pos_tx(2));
annotStrings{end+1} = sprintf('Est  TX: (%.2f, %.2f)', optCoord(1), optCoord(2));
annotStrings{end+1} = sprintf('Opt L: %.2f', L_peak);
annotStrings{end+1} = sprintf('Error: %.2f m', error);
annotStrings{end+1} = sprintf('SNR:  %.0f dB', SNR_dB(1,1));

%% Display the results
for i = 1:length(annotStrings)
    fprintf('%s\n', annotStrings{i});
end

%% Additional function to plot the 3D map


progressbar('reset', 1); % Reset progress bar
[X, Y, L] = l4c.calLikelihood4Area(area_size, pos_rx, rot_abs, y_centralised(1,:,:)', ELEMENT_NUM, nPower);
progressbar('end');  % This will display the total runtime

map3d.plots(X, Y, L, ...
    'Title', '3D Visualization of The Likelihood Function', ...
    'XLabel', 'x (m)', ...
    'YLabel', 'y (m)', ...
    'ZLabel', 'L(x,y)', ...
    'pos_tx', pos_tx, ...
    'pos_rx', pos_rx, ...
    'rot_abs', rot_abs, ...
    'area_size', area_size, ...
    'aoa_act', aoa_act, ...
    'angle_limit', ABS_ANGLE_LIM, ...
    'show_options', [SHOW_LIMITS, SHOW_EXTRA], ...
    'ShowAnnotation', true, ...
    'AnnotationStrings', annotStrings);
