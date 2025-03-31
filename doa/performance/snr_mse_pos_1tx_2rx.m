clear; clc; close all;
%% === User inputs
ITERATION = 1; % Number of Monte Carlo iterations
TIME_INST_NUM = 50; % Number of time instances
SNR_dB =repmat((0:1:20)', 1, 2); %dB
NUM_RX = 2; % Number of receivers
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree, which will 
RESOLUTION = 0.01; % Angle resolution in degreetranslate to -lim to lim (symmetric)
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 2; % Number of elements in the ULA
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
% --- Location Calculation
area_size = 100;   % 100x100 meter area
pos_tx = [50, 50]; % Transmitter position (x, y) in meters

n_param = length(SNR_dB); % Number of positions to test
progressbar('reset', n_param); % Reset progress bar
progressbar('displaymode', 'append'); % Reset progress bar
progressbar('minimalupdateinterval', 1); % Update progress bar every x seconds
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(NUM_RX, 1);  % W - Transmit signal power
sub_carrier = (1:NUM_RX)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA
%% === Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'ML_sync', 'MUSIC', 'BF', 'MVDR'} ... % extra args are defined later
);
num_methods = numel(doa_est_methods);  % Automatically get number of methods from struct array
% Preallocate MSE arrays
rmse_values = zeros(n_param, num_methods);
CRB_values = zeros(n_param, 1);
% square_err = zeros(ITERATION, num_methods);
map2d = Map2D();
%% === Loop through each SNR value
for snr_idx=1:n_param
    progressbar('step'); % Update progress bar
    % Pre-calculate required values outside loop
    tx_num = size(pos_tx, 1);
    % Generate base signal
    s_t = sqrt(P_t(NUM_RX)) .* exp(1j * 2 * pi * sub_carrier(NUM_RX) * t);
    %% ==== Define extra arguments for each method
    doa_est_methods(1).extra_args = {s_t};  % For ML_sync
    doa_est_methods(2).extra_args = {tx_num};  % For MUSIC
    doa_est_methods(3).extra_args = {};        % For MVDR
    doa_est_methods(4).extra_args = {};        % For BF
    % Calculate average energy of the signal
    if FIXED_TRANS_ENERGY == true
        % Average engery is fixed for whole transmission time
        % regardless of the number of time instances.
        avg_E = 1;
    else
        % Average engery is fixed for one time instance,
        % so the whole energy over all time instances need to be calculated.
        avg_E = avg_amp_gain^2 * P_t(NUM_RX) * T * Fs;
    end
    %% === Monte Carlo iterations
    for itr = 1:ITERATION
        % --- Location Refresh for each iteration
        aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], NUM_RX, 1); % true Angle of Arrival from RX to TX that will later be transformed to the absolute angle
        pos_rx = area_size*rand(NUM_RX, 2); % Random Receiver position (x, y) in meters
        angle_rx_tx_abs = zeros(NUM_RX, 1);
        for i = 1:NUM_RX
            % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
            angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
        end
        rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
        
        
        %% === DoA Estimation Algorithm
        aoa_rel_est = zeros(NUM_RX, num_methods);
        rays_abs = cell(NUM_RX, num_methods);
        % Initialize channel model
        channel = ChannelModels();
        for method_idx = 1:num_methods
            %% === Loop through each RX to find the ray from AoA
            for rx_idx=1:NUM_RX
                % Calculate noise parameters with the corresponding average energy and SNR
                nPower = avg_E/db2pow(SNR_dB(snr_idx, rx_idx));
                %% === Generate original signal received at Rx
                y_los = channel.LoS(s_t, avg_amp_gain);
                y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
                y_awgn = channel.AWGN(y_ula, nPower);
                %% === DoA Estimation Algorithm
                ula = ULA(lambda, ELEMENT_NUM, element_spacing);
                estimator = DoAEstimator(ula, sweeping_angle, aoa_act(rx_idx));
                result = estimator.(doa_est_methods(method_idx).name)(y_awgn, doa_est_methods(method_idx).extra_args{:});
                aoa_rel_est(rx_idx, method_idx) = result.aoa_est;
                rays_abs{rx_idx, method_idx} = map2d.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx, method_idx));
            end
            % --- Calculate the aoa intersection point and the RMSE for each method
            aoa_intersect = map2d.calDoAIntersect(rays_abs{1, method_idx}, rays_abs{2, method_idx});
            rmse_values(snr_idx, method_idx) = rmse_values(snr_idx, method_idx) + sqrt((pos_tx(1,1)-aoa_intersect.x)^2 + (pos_tx(1,2)-aoa_intersect.y)^2);
        end
    end
    rmse_values(snr_idx, :) = rmse_values(snr_idx, :) / ITERATION; % RMSE for current position
end

%% === Plotting
figure;
for i= 1:num_methods % Plot the the estimation methods' RMSE of the coordinates
    semilogy(mean(SNR_dB, 2), mean(rmse_values(:,i), 2), 'LineWidth', 1, 'DisplayName', strrep(doa_est_methods(i).name, '_', ' '));grid on; hold on;
end
title(['Coordinate RMSE of ', num2str(TIME_INST_NUM), 'ts']); legend("AutoUpdate","on");
xlabel('Signal to Noise Ratio (SNR)'); ylabel('Root Mean Square Error (RMSE)');
