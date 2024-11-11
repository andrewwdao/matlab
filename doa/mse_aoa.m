clear; clc; close all;
%% ====================== User inputs
ITERATION = 5;
TIME_INST_NUM = 10;
SNR_dB = 10; %dB

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
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_num = 2;   % Number of elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = (-90:1:90)'; % Angle range for sweeping to find the AoA

%% ==== Loop through each Tx position to test the accuracy from measuring the MSE
tic
mse_values = ones(size(tx_pos,1), 5);
for idx=1:size(tx_pos, 1)  % looping through the transmitter positions
    % Define the Channel and Antenna Array Model
    channel = ChannelModel(tx_pos{idx}, rx_pos, lambda, element_num, lambda/2);
    square_err = zeros(ITERATION, size(mse_values, 2)); % sync ML, Beamforming, MVDR, MUSIC
    for itr=1:ITERATION  % Iterating to get the average MSE
        %% ====================== Generate original signal received at Rx
        s_t = sqrt(P_t(idx)) .* exp(1j * 2 * pi * sub_carrier(idx) * t);
        y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
        y_ula = channel.applyULA(y_los);  % Apply ULA characteristics to the received signal
        % Calculate the energy of the whole signal transmitted to correct the noise power
        avg_E =  avg_amp_gain^2 * P_t(idx) * T * Fs; % average received signal energy, \sum_{t=1}^{T} ||s(t)||^2
        nPower = avg_E/db2pow(SNR_dB); % W - White noise power - noise variance
        y_awgn = channel.AWGN(y_ula, nPower);
        %% ---------------------- DoA Estimation Algorithm ----------------------
        estimator = DoAEstimator(y_awgn, size(tx_pos{idx},1), lambda, element_num, element_spacing, sweeping_angle);
        [est_aoa_sync, ~] = estimator.ML_sync(s_t);
        square_err(itr,1) = (est_aoa_sync - true_AoA(idx)).^2;
        [est_aoa_async, ~] = estimator.ML_async(s_t);
        square_err(itr,2) = (est_aoa_async - true_AoA(idx)).^2;
        [est_aoa_bf, ~] = estimator.BF();
        square_err(itr,3) = (est_aoa_bf - true_AoA(idx)).^2;
        [est_aoa_mvdr, ~] = estimator.MVDR();
        square_err(itr,4) = (est_aoa_mvdr - true_AoA(idx)).^2;
        [est_aoa_music, ~] = estimator.MUSIC();
        square_err(itr,5) = (est_aoa_music - true_AoA(idx)).^2;
    end
    for col = 1:size(mse_values, 2) % get the number of methods to mesure
        mse_values(idx, col) = mean(square_err(:, col));
    end
end
runtime = toc %#ok<NOPTS>
%% Plot the MSE
figure;
type = {'Sync ML', 'ASync ML', 'Beamforming', 'MVDR', 'MUSIC'};
for i= 1:size(mse_values, 2) % get the number of methods to mesure
    semilogy(true_AoA, mse_values(:,i), 'LineWidth', 1, 'DisplayName', type{i});
    grid on; hold on;
end
title(['Mean Square Error (MSE) for SNR=', num2str(SNR_dB), 'dB']); legend("AutoUpdate","on");
xlabel('Angle of Arrival (AoA)'); ylabel('MSE');
