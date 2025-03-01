clear; clc; close all;
%% === User inputs
ITERATION = 1; % Number of Monte Carlo iterations
TIME_INST_NUM = 1; % Number of time instances
TRUE_ANGLE = 0;
SNR_dB =(-20:2:20)'; %dB
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4;   % Number of elements in the ULA

%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
rx_pos = [10, 50;]; % Receiver position (x, y) in meters
pos_tx = [rx_pos(1) + 10*cosd(TRUE_ANGLE), rx_pos(2) + 10*sind(TRUE_ANGLE);]; % Transmitter position (x, y) in meters
aoa_act = zeros(size(rx_pos, 1), size(pos_tx, 1));
dist_act = zeros(size(rx_pos, 1), size(pos_tx, 1));
for i = 1:size(rx_pos, 1)
    for j = 1:size(pos_tx, 1)
        aoa_act(i,j) = atan2d(pos_tx(j,2) - rx_pos(i,2), pos_tx(j,1) - rx_pos(i,1)); % AoA in degrees - atan(y_tx-y_rx/x_tx-x_rx)
        dist_act = sqrt((pos_tx(j,1) - rx_pos(i,1))^2 + (pos_tx(j,2) - rx_pos(i,2))^2); % Euclidean distance between Tx and Rx - sqrt((x_tx-x_rx)^2 + (y_tx-y_rx)^2)
    end
end


n_param = length(SNR_dB); % Number of positions to test
progressbar('reset', n_param); % Reset progress bar
progressbar('displaymode', 'append'); % Reset progress bar
progressbar('minimalupdateinterval', 1); % Update progress bar every x seconds
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(size(pos_tx));  % W - Transmit signal power
sub_carrier = (1:size(pos_tx, 1))' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = (-90:0.1:90); % Angle range for finding the AoA
%% === Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'ML_sync', 'MUSIC', 'BF', 'MVDR'} ... % extra args are defined later
);
num_methods = numel(doa_est_methods);  % Automatically get number of methods from struct array
% Preallocate MSE arrays
mse_values = zeros(n_param, num_methods);
square_err = zeros(ITERATION, num_methods);
% Preallocate CRB array
CRB_values = zeros(n_param, 1);
CRB_Stoica_values = zeros(n_param, 1);
%% === Loop through each SNR value
for idx=1:n_param
    progressbar('advance'); % Update progress bar
    % Pre-calculate required values outside loop
    tx_num = size(pos_tx, 1);
    % Generate base signal
    s_t = sqrt(P_t(tx_num)) .* exp(1j * 2 * pi * sub_carrier(tx_num) * t);
    %% ==== Define extra arguments for each method
    doa_est_methods(1).extra_args = {s_t};  % For ML_sync
    doa_est_methods(2).extra_args = {tx_num};  % For MUSIC
    doa_est_methods(3).extra_args = {};        % For MVDR
    doa_est_methods(4).extra_args = {};        % For BF
    % Initialize channel model
    channel = ChannelModels();
    % Calculate average energy of the signal
    if FIXED_TRANS_ENERGY == true
        % Average engery is fixed for whole transmission time
        % regardless of the number of time instances.
        avg_E = 1;
    else
        % Average engery is fixed for one time instance,
        % so the whole energy over all time instances need to be calculated.
        avg_E = avg_amp_gain^2 * P_t(tx_num) * T * Fs;
    end
    % Calculate noise parameters with the corresponding average energy and SNR
    nPower = avg_E/db2pow(SNR_dB(idx));
    % Calculate CRB for the current position
    CRB_values(idx) = channel.CRB_det_1d_simp(s_t, nPower, aoa_act, ELEMENT_NUM, lambda);
    CRB_Stoica_values(idx) = channel.CRB_det_1d(s_t, nPower, aoa_act, ELEMENT_NUM, element_spacing, lambda);
    %% === Monte Carlo iterations
    for itr = 1:ITERATION
        % Generate received signal
        y_los = channel.LoS(s_t, avg_amp_gain);
        y_ula = channel.applyULA(y_los, aoa_act, ELEMENT_NUM, element_spacing, lambda);
        y_awgn = channel.AWGN(y_ula, nPower);
        %% === DoA Estimation Algorithm
        ula = ULA(lambda, ELEMENT_NUM, element_spacing);
        estimator = DoAEstimator(ula, sweeping_angle, aoa_act);
        for m = 1:num_methods
            result = estimator.(doa_est_methods(m).name)(y_awgn, doa_est_methods(m).extra_args{:});
            square_err(itr, m) = result.square_err;
        end
    end
    mse_values(idx, :) = mean(square_err, 1); % MSE for current position
end

%% === Plot the MSE vs SNR
figure;
semilogy(SNR_dB, CRB_values, 'r--', 'LineWidth', 1.5, 'DisplayName', 'CRB');grid on; hold on; % Plot CRB
semilogy(SNR_dB, CRB_Stoica_values, 'b--', 'LineWidth', 1.5, 'DisplayName', 'CRB Stoica');grid on; hold on; % Plot CRB
for i= 1:num_methods % Plot the rest of the estimation methods' MSE
    semilogy(SNR_dB, mse_values(:,i), 'LineWidth', 1, 'DisplayName', strrep(doa_est_methods(i).name, '_', ' '));
end
title(['MSE of ', num2str(TIME_INST_NUM), 'ts with total AoA=', num2str(TRUE_ANGLE), 'deg']); legend("AutoUpdate","on");
% ylim([1e-5, 1e4]);
xlabel('Signal to Noise Ratio (SNR)'); ylabel('MSE');
