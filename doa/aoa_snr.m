clear; clc; close all;
%% ====================== User inputs
ITERATION = 100;
TRUE_ANGLE = -30;
SNR_db =-50:10:50; %dB

%% ====================== Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
rx_pos = [10, 50;]; % Receiver position (x, y) in meters
tx_pos = [rx_pos(1) + 10*cosd(TRUE_ANGLE), rx_pos(2) + 10*sind(TRUE_ANGLE);]; % Transmitter position (x, y) in meters
element_num = 4;   % Number of elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:1:90; % Angle range for finding the AoA
P_t = 1; % W - Transmit signal power
%% ====================== Generate original signal received at Rx
s_t = sqrt(P_t) .* exp(1j * 2 * pi * 1000);  % Generate complex sinusoid transmitted signal
channel = ChannelModel(tx_pos, rx_pos, lambda, element_num, lambda/2);
array = ArrayModel(tx_pos, rx_pos, lambda, element_num, lambda/2);
y_los = channel.Direct(s_t);  % Received signal at the receiver
y_ula = array.applyULA(y_los);  % Apply ULA characteristics to the received signal
%% ====================== Loop through each SNR value
mse_values = zeros(length(SNR_db),2);
for idx=1:length(SNR_db)
    nPower = P_t/db2pow(SNR_db(idx)); % W - White noise power
    square_err = zeros(ITERATION, 2); % sync ML, Beamforming
    for itr=1:ITERATION  % Iterating to get the average MSE
        %% ---------------------- DoA Estimation Algorithm ----------------------
        % Add random noise to the received signal
        noise = sqrt(nPower/2) * (randn(size(y_ula)) + 1j * randn(size(y_ula))); % White noise
        y_t = y_ula + noise;
        % Start the estimation
        estimator = DoAEstimator(y_t, size(tx_pos,1), lambda, element_num, element_spacing, sweeping_angle);
        [est_aoa_sync, ~] = estimator.ML_sync();
        square_err(itr,1) = (est_aoa_sync - TRUE_ANGLE).^2;
        [est_aoa_bf, ~] = estimator.BF();
        square_err(itr,2) = (est_aoa_bf - TRUE_ANGLE).^2;
    end
    mse_values(idx,1) = mean(square_err(:,1));
    mse_values(idx,2) = mean(square_err(:,2));
end
figure; hold on; grid on;
plot(SNR_db, mse_values(:,1), 'r', 'LineWidth', 1, 'DisplayName', 'Sync ML');
plot(SNR_db, mse_values(:,2), 'g', 'LineWidth', 2, 'DisplayName', 'Beamforming');
title('Mean Square Error (MSE)'); legend("AutoUpdate","on");
xlabel('Signal to Noise Ratio (SNR)'); ylabel('MSE');
