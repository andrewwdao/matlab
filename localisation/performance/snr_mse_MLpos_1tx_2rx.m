clear; clc; close all;

%% User Inputs and Configurations
ITERATION = 500; % Number of Monte Carlo iterations
RX_NUM = 2;                         % Number of receivers
OPT_GRID_DENSITY = 10; % Define a coarse grid for initial guesses
SNR_dB = repmat((-10:3:20)', 1, 2);       % SNR in dB
ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 1;                   % Angle resolution (degrees)
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
pos_rx = [21, 51; 20, 40]; % 13
aoa_act = [0; 0];            % True AoA from Rx to Tx

n_param = length(SNR_dB); % Number of positions to test
progressbar('reset', ITERATION); % Reset progress bar

angle_rx_tx_abs = zeros(RX_NUM, 1);
for i = 1:RX_NUM
    % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
    angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
end
rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees

%% Signal and channel configurations
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(RX_NUM, 1);  % W - Transmit signal power
sub_carrier = (1:RX_NUM)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA

%% === Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'ML_sync'} ... % extra args are defined later
);
% doa_est_methods = struct(...
%     'name', {'ML_sync', 'MUSIC', 'MVDR', 'BF'} ... % extra args are defined later
% );
num_methods = numel(doa_est_methods);  % Automatically get number of methods from struct array
% Preallocate MSE arrays
rmse_values = zeros(n_param, num_methods);
rmse_values_ml = zeros(n_param, 1);
CRB_values = zeros(n_param, 1);
% square_err = zeros(ITERATION, num_methods);
% Initialize channel model
channel = ChannelModels();
map2d = Map2D();
l4c = Likelihood4Coordinates();
optimiser = gridOptimiser();
ula = ULA(lambda, ELEMENT_NUM, element_spacing);

%% Pre-calculate required values outside loop
tx_num = size(pos_tx, 1);
% Generate original transmitted signal
s_t = sqrt(P_t(RX_NUM)) .* exp(1j * 2 * pi * sub_carrier(RX_NUM) * t);
%% ==== Define extra arguments for each method
doa_est_methods(1).extra_args = {s_t};  % For ML_sync
% doa_est_methods(2).extra_args = {tx_num};  % For MUSIC
% doa_est_methods(3).extra_args = {};        % For MVDR
% doa_est_methods(4).extra_args = {};        % For BF
% Calculate average energy of the signal
avg_E = FIXED_TRANS_ENERGY * 1 + ~FIXED_TRANS_ENERGY * (avg_amp_gain^2 * P_t(RX_NUM) * T * Fs);
%% For each Rx, generate received signal
w = cell(RX_NUM, 1);
aoa_rel_est = zeros(RX_NUM, num_methods);
rays_abs = cell(RX_NUM, num_methods);
%% === Monte Carlo iterations
for itr = 1:ITERATION
    %% === Update progress bar
    progressbar('advance');
    %% --- Location Refresh for each iteration - ONLY ENABLE FOR RANDOMISED RX AND AOA
    aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1); % true Angle of Arrival from RX to TX that will later be transformed to the absolute angle
    pos_rx = area_size*rand(RX_NUM, 2); % Random Receiver position (x, y) in meters
    angle_rx_tx_abs = zeros(RX_NUM, 1);
    for i = 1:RX_NUM
        % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
        angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
    end
    rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
    %% === Loop through each SNR value
    for snr_idx=1:n_param
        %% === Loop through each method
        for method_idx = 1:num_methods
            %% === Loop through each RX to find the ray from AoA
            for rx_idx=1:RX_NUM
                % Calculate noise power with the corresponding average energy and SNR
                nPower = avg_E/db2pow(SNR_dB(snr_idx, rx_idx));
                %% === Generate signal received at Rx
                y_los = channel.LoS(s_t, avg_amp_gain);
                y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
                y_awgn = channel.AWGN(y_ula, nPower);
                %% === append received signal to a centralised array for direct ML estimation
                w{rx_idx} = y_awgn;
                %% === DoA Estimation Algorithm
                estimator = DoAEstimator(ula, sweeping_angle, aoa_act(rx_idx), 'opt', OPT_GRID_DENSITY);
                aoa_rel_est(rx_idx, method_idx) = estimator.(doa_est_methods(method_idx).name)(y_awgn, doa_est_methods(method_idx).extra_args{:}).aoa_est;
                rays_abs{rx_idx, method_idx} = map2d.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx, method_idx));
            end
            % --- Calculate the aoa intersection point and the RMSE for each method
            aoa_intersect = map2d.calDoAIntersect(rays_abs{1, method_idx}, rays_abs{2, method_idx});
            rmse_values(snr_idx, method_idx) = rmse_values(snr_idx, method_idx) + sqrt((pos_tx(1,1)-aoa_intersect.x)^2 + (pos_tx(1,2)-aoa_intersect.y)^2);
            %% === direct ML estimation
            objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
            [optCoord, ~] = optimiser.fmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
            rmse_values_ml(snr_idx, 1) = rmse_values_ml(snr_idx, 1) + sqrt((pos_tx(1,1)-optCoord(1))^2 + (pos_tx(1,2)-optCoord(2))^2);
        end
        rmse_values(snr_idx, :) = rmse_values(snr_idx, :) / ITERATION; % RMSE for DoA estimation methods
        rmse_values_ml(snr_idx, :) = rmse_values_ml(snr_idx, :) / (num_methods * ITERATION); % RMSE for ML, which is num_methods times more than the other methods
    end
end

%% === Plotting
figure;
for i= 1:num_methods % Plot the the estimation methods' RMSE of the coordinates
    semilogy(mean(SNR_dB, 2), mean(rmse_values(:,i), 2), 'LineWidth', 1, 'DisplayName', strrep(doa_est_methods(i).name, '_', ' '));grid on; hold on;
end
semilogy(mean(SNR_dB, 2), mean(rmse_values_ml(:,1), 2), 'LineWidth', 1, 'DisplayName', 'Maximum Likelihood');grid on; hold on;
title(['Coordinate RMSE of ', num2str(TIME_INST_NUM), 'ts']); legend("AutoUpdate","on");
xlabel('Signal to Noise Ratio (SNR)'); ylabel('Root Mean Square Error (RMSE)');
