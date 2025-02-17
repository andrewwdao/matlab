clear; clc; close all;
%% === User inputs
ITERATION = 100; % Number of Monte Carlo iterations
SNR_dB = -30; %dB
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree
TIME_INST_NUM = 150; % Number of time instances
aoa_act = 30; % True Angle of Arrival
RESOLUTION = 0.1; % Angle resolution in degree
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4; % Number of elements in the ULA
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
pos_rx = [0, 50;]; % Receiver position (x, y) in meters
pos_tx = [pos_rx(1)+40*cosd(aoa_act), pos_rx(2)+40*sind(aoa_act);]; % Transmitter position (x, y) in meters

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
    avg_E = avg_amp_gain^2 * P_t(RX_NUM) * T * Fs;
end
% Calculate noise parameters with the corresponding average energy and SNR
nPower = avg_E/db2pow(SNR_dB);
progressbar('reset', ITERATION); % Reset progress bar
progressbar('displaymode', 'append'); % Reset progress bar
progressbar('minimalupdateinterval', 1); % Update progress bar every x seconds
est_aoas = zeros(ITERATION, 1);
parfor itr=1:ITERATION
    progressbar('advance'); % Update progress bar
    % Initialize channel model
    channel = ChannelModels();
    y_los = channel.LoS(s_t, avg_amp_gain);
    y_ula = channel.applyULA(y_los, aoa_act, ELEMENT_NUM, element_spacing, lambda);
    y_awgn = channel.AWGN(y_ula, nPower);

    %% === DoA Estimation Algorithm
    estimator_angle = DoAEstimator(y_awgn, size(pos_tx,1), lambda, ELEMENT_NUM, element_spacing, sweeping_angle, aoa_act);
    est_aoas(itr) = estimator_angle.ML_sync(s_t).aoa_est;
    % result = estimator_angle.BF();
    % result = estimator_angle.MVDR();
    % result = estimator_angle.MUSIC();
end

%% === Calculate statistics
mean_est_aoa = mean(est_aoas);
std_est_aoa = std(est_aoas);

%% === Plotting the Empirical PDF of Estimated Angle
figure; hold on; grid on;
histogram(est_aoas, 'Normalization', 'pdf');
xlabel('Estimated Angle (degrees)');
ylabel('Probability Density');
title(sprintf('Empirical PDF of AoA for %d iterations\nMean: %.2f°, Std: %.2f°', ITERATION, mean_est_aoa, std_est_aoa));
xlim([20 40]);
% Add vertical line for true angle
xline(aoa_act, 'r--', 'True Angle', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');