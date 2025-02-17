clear; clc; close all;
%% === User inputs
ITERATION = 1; % Number of Monte Carlo iterations
TIME_INST_NUM = 50; % Number of time instances
SNR_dB =repmat((0:1:20)', 1, 2); %dB
RX_NUM = 2; % Number of receivers
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
    'name', {'BF', 'MVDR', 'MUSIC'}, ...
    'transmitted_signal_required', {false, false, false});
    % 'name', {'ML_sync', 'BF'}, ..., 'MVDR', 'MUSIC'}, ...
    % 'transmitted_signal_required', {true, false});%, false, false});
num_methods = numel(doa_est_methods);  % Automatically get number of methods from struct array
% Preallocate MSE arrays
rmse_values = zeros(n_param, num_methods);
CRB_values = zeros(n_param, 1);
% square_err = zeros(ITERATION, num_methods);
estimator_coor = PosEstimator2D();
%% === Loop through each SNR value
for snr_idx=1:n_param
    progressbar('advance'); % Update progress bar
    % Pre-calculate required values outside loop
    % Generate base signal
    s_t = sqrt(P_t(RX_NUM)) .* exp(1j * 2 * pi * sub_carrier(RX_NUM) * t);
    % Calculate average energy of the signal
    if FIXED_TRANS_ENERGY == true
        % Average engery is fixed for whole transmission time
        % regardless of the number of time instances.
        avg_E = 1;
    else
        % Average engery is fixed for one time instance,
        % so the whole energy over all time instances need to be calculated.
        avg_E = avg_amp_gain^2 * P_t(RX_NUM) * T * Fs;
    end
    %% === Monte Carlo iterations
    for itr = 1:ITERATION
        % --- Location Refresh for each iteration
        aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1); % true Angle of Arrival from RX to TX that will later be transformed to the absolute angle
        pos_rx = area_size*rand(RX_NUM, 2); % Random Receiver position (x, y) in meters
        angle_rx_tx_abs = zeros(RX_NUM, 1);
        for i = 1:RX_NUM
            % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
            angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
        end
        rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
        
        
        %% === DoA Estimation Algorithm
        aoa_rel_est = zeros(RX_NUM, num_methods);
        rays_abs = cell(RX_NUM, num_methods);
        % Initialize channel model
        channel = ChannelModels();
        for method_idx = 1:num_methods
            %% === Loop through each RX to find the ray from AoA
            for rx_idx=1:RX_NUM
                % Calculate noise parameters with the corresponding average energy and SNR
                nPower = avg_E/db2pow(SNR_dB(snr_idx, rx_idx));
                %% === Generate original signal received at Rx
                y_los = channel.LoS(s_t, avg_amp_gain);
                y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
                y_awgn = channel.AWGN(y_ula, nPower);

                estimator_angle = DoAEstimator(y_awgn, size(pos_tx,1), lambda, ...
                    ELEMENT_NUM, element_spacing, sweeping_angle, aoa_act(rx_idx));
                if doa_est_methods(method_idx).transmitted_signal_required
                    aoa_rel_est(rx_idx, method_idx) = estimator_angle.(doa_est_methods(method_idx).name)(s_t).aoa_est;
                else
                    aoa_rel_est(rx_idx, method_idx) = estimator_angle.(doa_est_methods(method_idx).name)().aoa_est;
                end
                rays_abs{rx_idx, method_idx} = estimator_coor.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx, method_idx));
            end
            % --- Calculate the aoa intersection point and the RMSE for each method
            aoa_intersect = estimator_coor.calDoAIntersect(rays_abs{1, method_idx}, rays_abs{2, method_idx});
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
