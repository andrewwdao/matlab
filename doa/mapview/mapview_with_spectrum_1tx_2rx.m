clear; clc; close all;
%% === User inputs
SNR_dB = 20; %dB
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree
TIME_INST_NUM = 100; % Number of time instances
aoa_act = [30, 35]; % True Angle of Arrival
RESOLUTION = 0.1; % Angle resolution in degree
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 16; % Number of elements in the ULA

SHOW_LIMITS = true; % Show the detecting limits of the RXs (with known limitation)
SHOW_EXTRA = false; % Show extra information such as the AoA and the intersection point
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
area_size = 100;   % 100x100 meter area
pos_tx = zeros(length(aoa_act), 2);
for i = 1:length(aoa_act)
    pos_tx(i, :) = [0+40*cosd(aoa_act(i)), 50+40*sind(aoa_act(i));]; % Transmitter position (x, y) in meters
end
% pos_tx = [0+40*cosd(aoa_act), 50+40*sind(aoa_act);]; % Transmitter position (x, y) in meters
pos_rx = [0, 50;]; % Receiver position (x, y) in meters
tx_num = size(pos_tx, 1);
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


%% === DoA Estimation Algorithm
% Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'MUSIC', 'MVDR', 'BF'}, ...
    'extra_args', {{tx_num}, {}, {}} ...
);
method_num = numel(doa_est_methods);  % Automatically get number of methods from struct array
spectrum_cell = cell(method_num, tx_num);
aoa_est_cell = cell(method_num, tx_num);
map2d = Map2D();
channel = ChannelModels();  % Initialize channel model
for method_idx=1:method_num
    for tx_idx = 1:tx_num
        y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
        y_ula = channel.applyULA(y_los, aoa_act(tx_idx), ELEMENT_NUM, element_spacing, lambda);  % Apply ULA characteristics to the received signal
        y_awgn = channel.AWGN(y_ula, nPower);
        %% === DoA Estimation Algorithm
        ula = ULA(lambda, ELEMENT_NUM, element_spacing);
        estimator = DoAEstimator(ula, sweeping_angle, aoa_act(tx_idx));
        result = estimator.(doa_est_methods(method_idx).name)(y_awgn, doa_est_methods(method_idx).extra_args{:});
        spectrum_cell{method_idx, tx_idx} = result.spectrum_dB;
        aoa_est_cell{method_idx, tx_idx} = result.aoa_est;
    end
end

%% === Plotting
map2d.plotDetailed(pos_tx, pos_rx, 0, area_size, aoa_act, ABS_ANGLE_LIM, [SHOW_LIMITS, SHOW_EXTRA], sweeping_angle, spectrum_cell, {doa_est_methods.name}, aoa_est_cell);
