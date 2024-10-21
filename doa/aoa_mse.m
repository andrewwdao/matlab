clear; clc; close all;
%% ====================== User inputs
ITERATION = 100;
SNR_db =-50; %dB

%% ====================== Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
true_AoA = (-90:10:90)'; % Angle range for sweeping to find the AoA
area_size = 100;   % 100x100 meter area
rx_pos = [10,50;]; % Receiver position (x, y) in meters
tx_pos = cell(size(true_AoA));
for i = 1:length(true_AoA)
    tx_pos{i} = [rx_pos(1) + 10*cosd(true_AoA(i)), rx_pos(2) + 10*sind(true_AoA(i))];
end
element_num = 128;   % Number of elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = (-90:2:90)'; % Angle range for sweeping to find the AoA
P_t = ones(size(tx_pos));    % W - Transmit signal power
tic
%% ====================== Loop through each Tx position to test the accuracy from measuring the MSE
mse_values = ones(size(tx_pos,1),4);
for idx=1:size(tx_pos, 1)  % looping through the transmitter positions
    % Define the Channel and Antenna Array Model
    channel = ChannelModel(tx_pos{idx}, rx_pos, lambda, element_num, lambda/2);
    square_err = zeros(ITERATION, 4); % sync ML, Beamforming, MVDR, MUSIC
    for itr=1:ITERATION  % Iterating to get the average MSE
        %% ====================== Generate original signal received at Rx
        s_t = sqrt(P_t(idx, 1)) .* exp(1j * 2 * pi * 1000);  % Generate complex sinusoid transmitted signal
        y_los = channel.LoS(s_t, 1);  % Received signal at the receiver
        y_ula = channel.applyULA(y_los);  % Apply ULA characteristics to the received signal
        % Calculate the noise power and the corresponding noise
        nPower = 1/db2pow(SNR_db); % W - White noise power
        y_awgn = channel.AWGN(y_ula, nPower);
        %% ---------------------- DoA Estimation Algorithm ----------------------
        estimator = DoAEstimator(y_awgn, size(tx_pos{idx},1), lambda, element_num, element_spacing, sweeping_angle);
        [est_aoa_sync, ~] = estimator.ML_sync();
        square_err(itr,1) = (est_aoa_sync - true_AoA(idx)).^2;
        [est_aoa_bf, ~] = estimator.BF();
        square_err(itr,2) = (est_aoa_bf - true_AoA(idx)).^2;
        % [est_aoa_mvdr, ~] = estimator.MVDR();
        % square_err(itr,3) = (est_aoa_mvdr - sweeping_angle(idx)).^2;
        % [est_aoa_music, ~] = estimator.MUSIC();
        % square_err(itr,4) = (est_aoa_music - sweeping_angle(idx)).^2;
    end
    mse_values(idx,1) = mean(square_err(:,1));
    mse_values(idx,2) = mean(square_err(:,2));
    % mse_values(idx,3) = mean(square_err(:,3));
    % mse_values(idx,4) = mean(square_err(:,4));
end
runtime = toc
%% Plot the MSE
figure; 
semilogy(true_AoA, mse_values(:,1), 'r', 'LineWidth', 1.25, 'DisplayName', 'Sync ML');
grid on; hold on;
semilogy(true_AoA, mse_values(:,2), 'g', 'LineWidth', 1.75, 'DisplayName', 'Beamforming');
% semilogy(sweeping_angle, mse_values(:,3), 'b', 'LineWidth', 0.5, 'DisplayName', 'MVDR');
% semilogy(sweeping_angle, mse_values(:,4), 'm', 'LineWidth', 0.25, 'DisplayName', 'MUSIC');
title('Mean Square Error (MSE)'); legend("AutoUpdate","on");
xlabel('Angle of Arrival (AoA)'); ylabel('MSE');

% %% --- Plotting
% vis_music = DoAVisualisation(...
%     "MUSIC", tx_pos, rx_pos, area_size, ...
%     sweeping_angle, spectrum_dB_music, ...
%         est_aoa_music);
% vis_music.plot();
% vis_mvdr = DoAVisualisation(...
%     "MVDR", tx_pos, rx_pos, area_size, ...
%     sweeping_angle, spectrum_dB_mvdr, ...
%         est_aoa_mvdr);
% vis_mvdr.plot();
% vis_beam = DoAVisualisation(...
%     "Beamforming", tx_pos, rx_pos, area_size, ...
%     sweeping_angle, spectrum_dB_beam, ...
%         est_aoa_beam);
% vis_beam.plot();
% vis_simple = DoAVisualisation(...
%     "Beamforming", tx_pos, rx_pos, area_size, ...
%     sweeping_angle, spectrum_dB_simple, ...
%         est_aoa_simple);
% vis_simple.plot();
