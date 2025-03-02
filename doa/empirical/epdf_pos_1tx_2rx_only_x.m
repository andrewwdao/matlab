clear; clc; close all;
%% === User inputs
ITERATION = 100; % Number of Monte Carlo iterations
SNR_dB = -120; %dB
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree
TIME_INST_NUM = 150; % Number of time instances
TRUE_ANGLE = 30; % True Angle of Arrival
RESOLUTION = 0.1; % Angle resolution in degree
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4; % Number of elements in the ULA
OPT_GRID_DENSITY = 3;

%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area

pos_tx = [50, 50]; % Transmitter position (x, y) in meters
pos_rx = [15.98, 54.27; 67.64, 69.28];
aoa_act = [-6.8; 45.6];

RX_NUM = size(pos_rx, 1); % Number of receivers
angle_rx_tx_abs = zeros(RX_NUM, 1);
for i = 1:RX_NUM  % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
    angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
end
rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees

avg_amp_gain = 1; % Average gain of the channel
P_t = 1;  % W - Transmit signal power1000Hz
sub_carrier = (1:1)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA
%% === Generate original signal received at Rx
s_t = sqrt(P_t(1)) .* exp(1j * 2 * pi * sub_carrier(1) * t);
% Calculate average energy of the signal
if FIXED_TRANS_ENERGY == true
    % Average engery is fixed for whole transmission time
    % regardless of the number of time instances.
    avg_E = 1;
else
    % Average engery is fixed for one time instance,
    % so the whole energy over all time instances need to be calculated.
    avg_E = avg_amp_gain^2 * P_t(1) * T * Fs;
end
% Calculate noise parameters with the corresponding average energy and SNR
nPower = avg_E/db2pow(SNR_dB);
progressbar('reset', ITERATION); % Reset progress bar
progressbar('displaymode', 'append'); % Reset progress bar
progressbar('minimalupdateinterval', 1); % Update progress bar every x seconds
rmse_values = zeros(ITERATION, 1);
map2d = Map2D();
% Initialize channel model
channel = ChannelModels();
for itr=1:ITERATION
    progressbar('advance'); % Update progress bar
    %% === Loop through each RX
    aoa_rel_est = zeros(RX_NUM, 1);
    rays_abs = cell(RX_NUM, 1);
    for rx_idx=1:RX_NUM
        y_los = channel.LoS(s_t, avg_amp_gain);
        y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
        y_awgn = channel.AWGN(y_ula, nPower);

        %% === DoA Estimation Algorithm
        ula = ULA(lambda, ELEMENT_NUM, element_spacing);
        estimator = DoAEstimator(ula, sweeping_angle, aoa_act(rx_idx), 'opt', OPT_GRID_DENSITY);
        result = estimator.ML_sync(y_awgn, s_t);
        % result = estimator.BF(y_awgn);
        % result = estimator.MVDR(y_awgn);
        % result = estimator.MUSIC(y_awgn, length(pos_tx));
        aoa_rel_est(rx_idx, 1) = result.aoa_est;
        rays_abs{rx_idx, 1} = map2d.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx, 1));
    end
    % --- Calculate the aoa intersection point and the RMSE for each method
    aoa_intersect = map2d.calDoAIntersect(rays_abs{1, 1}, rays_abs{2, 1});
    rmse_values(itr, 1) = round((pos_tx(1,1)-aoa_intersect.x), 0);
end

%% === Calculate statistics
mean_est_rmse = mean(rmse_values);
std_est_rmse = std(rmse_values);

%% === Plotting the Empirical PDF of Estimated Angle
figure; hold on; grid on;
histogram(rmse_values, 'BinWidth', 10, ...
    'Normalization', 'pdf', ... % Normalize to probability density
    'FaceColor', 'b', ...
    'EdgeColor', 'k');
xlabel('Estimated RMSE for Coordinate (m)');
ylabel('Probability Density');
title(sprintf('Empirical PDF of |x-x^| for %d iterations at %d dB SNR\nMean: %.2f m, Std: %.2f m', ITERATION, SNR_dB, mean_est_rmse, std_est_rmse));
xlim([-500 500]);
% Add vertical line for true angle
xline(0, 'r--', 'Optimal value', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
