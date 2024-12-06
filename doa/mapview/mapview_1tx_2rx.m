clear; clc; close all;
%% === User inputs
RX_NUM = 2; % Number of receivers
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree, which will 
TIME_INST_NUM = 1000; % Number of time instances
RESOLUTION = 0.1; % Angle resolution in degreetranslate to -lim to lim (symmetric)
SNR_dB = 100*ones(RX_NUM, 1); %dB
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4; % Number of elements in the ULA
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
% --- Location Calculation
area_size = 100;   % 100x100 meter area
true_aoa = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1); % true Angle of Arrival from RX to TX that will later be transformed to the absolute angle
tx_pos = [50, 50]; % Transmitter position (x, y) in meters
rx_pos = area_size*rand(RX_NUM, 2); % Random Receiver position (x, y) in meters
abs_angles = zeros(RX_NUM, 1);
for i = 1:RX_NUM
    % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
    abs_angles(i) = atan2d(tx_pos(2)-rx_pos(i,2), tx_pos(1)-rx_pos(i,1));
end
abs_rot = abs_angles - true_aoa; % Absolute rotation of the receiver in degrees

progressbar('reset', RX_NUM); % Reset progress bar
progressbar('displaymode', 'append'); % Reset progress bar
progressbar('minimalupdateinterval', 1); % Update progress bar every x seconds
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(RX_NUM, 1);  % W - Transmit signal power
sub_carrier = (1:RX_NUM)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA
%% === Loop through each RX
est_rel_angle = zeros(RX_NUM, 1);
abs_ray = cell(RX_NUM, 1);
coor_estimator = PosEstimator2D();
for rx_idx=1:RX_NUM
    progressbar('advance'); % Update progress bar
    % Pre-calculate required values outside loop
    % Generate base signal
    s_t = sqrt(P_t(RX_NUM)) .* exp(1j * 2 * pi * sub_carrier(RX_NUM) * t);
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
    nPower = avg_E/db2pow(SNR_dB(rx_idx));
    % Initialize channel model
    channel = ChannelModelAoA(true_aoa(rx_idx), lambda, ELEMENT_NUM, element_spacing);
    %% === Generate original signal received at Rx
    y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
    y_ula = channel.applyULA(y_los);  % Apply ULA characteristics to the received signal
    y_awgn = channel.AWGN(y_ula, nPower);

    %% === DoA Estimation Algorithm
    angle_estimator = DoAEstimator(y_awgn, size(tx_pos,1), lambda, ELEMENT_NUM, element_spacing, sweeping_angle, true_aoa);
    est_rel_angle(rx_idx) = angle_estimator.ML_sync(s_t).est_aoa;
    abs_ray{rx_idx} = coor_estimator.calAbsRay(rx_pos(rx_idx,:), tx_pos, abs_rot(rx_idx), est_rel_angle(rx_idx));
end
intersection = coor_estimator.calIntersection(abs_ray{1}, abs_ray{2});
%% === Plotting
figure('Name', 'Map Visualisation', 'WindowState', 'maximized'); clf; hold on;
vis_sync = MapVisual("ML Sync", tx_pos, rx_pos, area_size, ...
    sweeping_angle, [], est_rel_angle);
vis_sync.plotMapOnly(area_size, tx_pos, rx_pos, abs_ray, intersection);
