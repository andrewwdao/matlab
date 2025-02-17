clear; clc; close all;
%% === User inputs
RX_NUM = 2; % Number of receivers
SNR_dB = 100*ones(RX_NUM, 1); %dB
SHOW_LIMITS = true; % Show the detecting limits of the RXs (with known limitation)
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree
TIME_INST_NUM = 150; % Number of time instances
RESOLUTION = 0.1; % Angle resolution in degreetranslate to -lim to lim (symmetric)
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4; % Number of elements in the ULA
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
% --- Location Calculation
area_size = 100;   % 100x100 meter area
aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1); % true Angle of Arrival from RX to TX that will later be transformed to the absolute angle
% aoa_act = [-6.8; 45.6];
pos_tx = [50, 50]; % Transmitter position (x, y) in meters
pos_rx = area_size*rand(RX_NUM, 2); % Random Receiver position (x, y) in meters
% pos_rx = [15.98, 54.27; 67.64, 69.28];
angle_rx_tx_abs = zeros(RX_NUM, 1);
for i = 1:RX_NUM  % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
    angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
end
rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees

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
aoa_rel_est = zeros(RX_NUM, 1);
rays_abs = cell(RX_NUM, 1);
estimator_coor = PosEstimator2D();
% Initialize channel model
channel = ChannelModels();
for rx_idx=1:RX_NUM
    progressbar('advance'); % Update progress bar
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
    %% === Generate original signal received at Rx
    y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
    y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);  % Apply ULA characteristics to the received signal
    y_awgn = channel.AWGN(y_ula, nPower);

    %% === DoA Estimation Algorithm
    estimator_angle = DoAEstimator(y_awgn, size(pos_tx,1), lambda, ELEMENT_NUM, element_spacing, sweeping_angle, aoa_act);
    aoa_rel_est(rx_idx) = estimator_angle.ML_sync(s_t).aoa_est;
    rays_abs{rx_idx} = estimator_coor.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx), ABS_ANGLE_LIM);
end
aoa_intersect = estimator_coor.calDoAIntersect(rays_abs{1}, rays_abs{2});
%% === Plotting
fprintf('SNR (dB):\n')
fprintf('%.0f  ', SNR_dB);
fprintf('\n');
fprintf('RX Positions:\n');
for idx = 1:size(pos_rx, 1)
    fprintf('  Rx %d: x = %.2f, y = %.2f\n', idx, pos_rx(idx, 1), pos_rx(idx, 2));
end

fprintf('Absolute RX Rotations (degrees):\n  ');
fprintf('%.2f  ', rot_abs);
fprintf('\n');

fprintf('Absolute RX-TX Angles (degrees):\n  ');
fprintf('%.2f  ', angle_rx_tx_abs);
fprintf('\n');

fprintf('True AoA (degrees):\n  ');
fprintf('%.2f  ', aoa_act);
fprintf('\n');

fprintf('Estimated AoA (degrees):\n  ');
fprintf('%.2f  ', aoa_rel_est);
fprintf('\n');

fprintf('DoA Intersection Point (x, y):\n  %.2f, %.2f\n', aoa_intersect.x, aoa_intersect.y);
figure('Name', 'Map Visualisation'); clf; hold on;%, 'WindowState', 'maximized'); clf; hold on;
vis_sync = MapVisual("ML Sync", pos_tx, pos_rx, area_size, ...
    sweeping_angle, [], aoa_rel_est);
vis_sync.plotMapOnly(area_size, pos_tx, pos_rx, rays_abs, aoa_intersect, SHOW_LIMITS);
