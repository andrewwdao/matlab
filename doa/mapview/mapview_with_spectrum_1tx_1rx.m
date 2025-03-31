clear; clc; close all;
%% === User inputs
SNR_dB = 10; %dB
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree
TIME_INST_NUM = 150; % Number of time instances
aoa_act = 30; % True Angle of Arrival
RESOLUTION = 0.1; % Angle resolution in degree
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4; % Number of elements in the ULA


SHOW_LIMITS = true; % Show the detecting limits of the RXs (with known limitation)
SHOW_EXTRA = true; % Show extra information such as the AoA and the intersection point
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
pos_tx = [0+40*cosd(aoa_act), 50+40*sind(aoa_act);]; % Transmitter position (x, y) in meters
pos_rx = [0, 50;]; % Receiver position (x, y) in meters

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
    avg_E = avg_amp_gain^2 * P_t(NUM_RX) * T * Fs;
end
% Calculate noise parameters with the corresponding average energy and SNR
nPower = avg_E/db2pow(SNR_dB);

% Initialize channel model
channel = ChannelModels();
y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
y_ula = channel.applyULA(y_los, aoa_act, ELEMENT_NUM, element_spacing, lambda);
y_awgn = channel.AWGN(y_ula, nPower);

%% === DoA Estimation Algorithm
ula = ULA(lambda, ELEMENT_NUM, element_spacing);
estimator = DoAEstimator(ula, sweeping_angle, aoa_act);
% result = estimator.ML_sync(y_awgn, s_t);
result = estimator.BF(y_awgn);
% result = estimator.MVDR(y_awgn);
% result = estimator.MUSIC(y_awgn, length(pos_tx));

%% === Plotting
% figure; hold on; grid on;
% plot(result_sync.spectrum_dB, 'LineWidth', 1, 'DisplayName', "Sync");
% plot(result_bf.spectrum_dB, 'LineWidth', 1, 'DisplayName', "BF");
% title("Spectrum"); legend("AutoUpdate","on");
map2d = Map2D();
map2d.plotDetailed(pos_tx, pos_rx, 0, area_size, aoa_act, ABS_ANGLE_LIM, [SHOW_LIMITS, SHOW_EXTRA], sweeping_angle, {result.spectrum_dB}, "ML sync", {result.aoa_est});
