clear; clc; close all;
%% === User inputs
SNR_dB = 100; %dB
SHOW_LIMITS = true; % Show the detecting limits of the RXs (with known limitation)
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree
TIME_INST_NUM = 1000; % Number of time instances
TRUE_ANGLE = 30; % True Angle of Arrival
RESOLUTION = 0.1; % Angle resolution in degree
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4; % Number of elements in the ULA
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
tx_pos = [0+40*cosd(TRUE_ANGLE), 50+40*sind(TRUE_ANGLE);]; % Transmitter position (x, y) in meters
rx_pos = [0, 50;]; % Receiver position (x, y) in meters

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
% Initialize channel model
channel = ChannelModel(tx_pos, rx_pos, lambda, ELEMENT_NUM, lambda/2);
y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
y_ula = channel.applyULA(y_los);  % Apply ULA characteristics to the received signal
% Calculate the energy of the whole signal transmitted to correct the noise power
y_awgn = channel.AWGN(y_ula, nPower);

%% === DoA Estimation Algorithm
estimator = DoAEstimator(y_awgn, size(tx_pos,1), lambda, ELEMENT_NUM, element_spacing, sweeping_angle, TRUE_ANGLE);
result = estimator.ML_sync(s_t);
% result = estimator.BF();
% result = estimator.MVDR();
% result = estimator.MUSIC();
coor_estimator = PosEstimator2D();
abs_ray = cell(1, 1);
abs_ray{1} = coor_estimator.calAbsRays(rx_pos, tx_pos, 0, result.est_aoa, ABS_ANGLE_LIM);

%% === Plotting
% figure; hold on; grid on;
% plot(result_sync.spectrum_dB, 'LineWidth', 1, 'DisplayName', "Sync");
% plot(result_bf.spectrum_dB, 'LineWidth', 1, 'DisplayName', "BF");
% title("Spectrum"); legend("AutoUpdate","on");
vis_sync = MapVisual(...
    "ML Sync", tx_pos, rx_pos, area_size, ...
    sweeping_angle, result.spectrum_dB, ...
        result.est_aoa);
vis_sync.plotSingle(abs_ray, SHOW_LIMITS);

