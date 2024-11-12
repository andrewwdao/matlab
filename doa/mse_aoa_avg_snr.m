clear; clc; close all;
%% ====================== User inputs
ITERATION = 50;
TIME_INST_num = 15;
SNR_dB = 10; %dB - average SNR for all time instances
element_num = 4;   % Number of elements in the ULA

%% ====================== Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
true_AoA = (-85:5:85)'; % Angle range for sweeping to find the AoA
rx_pos = [10,50;]; % Receiver position (x, y) in meters
tx_pos = cell(size(true_AoA));
for i = 1:length(true_AoA)
    tx_pos{i} = [rx_pos(1) + 10*cosd(true_AoA(i)), rx_pos(2) + 10*sind(true_AoA(i))];
end
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(size(tx_pos));  % W - Transmit signal power
sub_carrier = (1:size(tx_pos,1))' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_num/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = (-90:1:90)'; % Angle range for sweeping to find the AoA

%% ==== Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'ML_sync', 'ML_async', 'BF', 'MVDR', 'MUSIC'}, ...
    'transmitted_signal_required', {true, true, false, false, false});
% Preallocate arrays
num_methods = numel(doa_est_methods);  % Automatically get number of methods from struct array
mse_values = zeros(size(tx_pos,1), num_methods);
square_err = zeros(ITERATION, num_methods);

%% ==== Loop through each Tx position to test the accuracy from measuring the MSE
tic; % Start timing the simulation
for idx = 1:size(tx_pos, 1)
    % Pre-calculate required values outside parfor loop
    tx_num = size(tx_pos{idx}, 1);
    % Generate base signal
    s_t = sqrt(P_t(idx)) .* exp(1j * 2 * pi * sub_carrier(idx) * t);
    % Initialize channel model
    channel = ChannelModel(tx_pos{idx}, rx_pos, lambda, element_num, element_spacing);
    % Pre-calculate noise parameters
    avg_E = avg_amp_gain^2 * P_t(idx) * T * Fs;
    nPower = avg_E/db2pow(SNR_dB);
    
    % Monte Carlo iterations
    for itr = 1:ITERATION
        % Generate received signal
        y_los = channel.LoS(s_t, avg_amp_gain);
        y_ula = channel.applyULA(y_los);
        y_awgn = channel.AWGN(y_ula, nPower);
        % Create local estimator for this iteration
        estimator = DoAEstimator(y_awgn, tx_num, lambda, ...
            element_num, element_spacing, sweeping_angle, channel.act_aoa);
        
        for m = 1:num_methods
            if doa_est_methods(m).transmitted_signal_required
                result = estimator.(doa_est_methods(m).name)(s_t);
            else
                result = estimator.(doa_est_methods(m).name)();
            end
            square_err(itr, m) = result.square_err;
        end
    end
    mse_values(idx, :) = mean(square_err, 1); % MSE for current position
end
runtime = toc %#ok<NOPTS>
%% Plot the MSE
figure;
for i= 1:num_methods % get the number of methods to mesure
    semilogy(true_AoA, mse_values(:,i), 'LineWidth', 1, 'DisplayName', strrep(doa_est_methods(i).name, '_', ' '));
    grid on; hold on;
end
title(['MSE of ', num2str(TIME_INST_num), 'ts with avg. SNR=', num2str(SNR_dB), 'dB']); legend("AutoUpdate","on");
xlabel('Angle of Arrival (AoA)'); ylabel('Mean Square Error (MSE)');
