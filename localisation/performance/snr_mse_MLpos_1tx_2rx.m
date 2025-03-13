clear; clc; close all;

%% User Inputs and Configurations
ITERATION = 5000; % Number of Monte Carlo iterations
OPT_GRID_DENSITY = 10; % Define a coarse grid for initial guesses
ABS_ANGLE_LIM = 1;                 % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 0.1;                   % Angle resolution (degrees)
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements
DOA_MODE = 'sweep';                % DoA estimation mode ('sweep' or 'opt')
SAFETY_DISTANCE = 2;             % Minimum distance between TX and RX (meters)
% Physical constants and wavelength
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength
% Transmitter, receiver positions angles
area_size = 100;
pos_tx = [50, 50];                  % Tx at center

RX_NUM = 2;                         % Number of receivers
SNR_dB = repmat((-10:2:20)', 1, RX_NUM);       % SNR in dB
% --- Fixed rx and aoa - optimal scenario
pos_rx = [21, 51; 20, 40]; % 13
aoa_act = [0; 0];            % True AoA from Rx to Tx

% --- Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
angle_rx_tx_abs = zeros(RX_NUM, 1);
for i = 1:RX_NUM
    % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
    angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
end
rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees

n_param = length(SNR_dB); % Number of positions to test
fprintf('Running Monte Carlo simulation with %d iterations...\n', ITERATION);
progressbar('reset', ITERATION*n_param); % Reset progress bar

%% Signal and channel configurations
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(RX_NUM, 1);  % W - Transmit signal power
sub_carrier = (1:RX_NUM)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA

%% Pre-calculate required values outside loop
tx_num = size(pos_tx, 1);
% Generate original transmitted signal
s_t = sqrt(P_t(RX_NUM)) .* exp(1j * 2 * pi * sub_carrier(RX_NUM) * t);
% Calculate average energy of the signal
avg_E = FIXED_TRANS_ENERGY * 1 + ~FIXED_TRANS_ENERGY * (avg_amp_gain^2 * P_t(RX_NUM) * T * Fs);
%% === Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'ML_async'}, ... % estimator methods
    'extra_args', {{}} ...% extra args required for specific type of estimator
);
num_methods = numel(doa_est_methods);  % Automatically get number of methods from struct array

%% Initialise classes and arrays
channel = ChannelModels();
y_los = channel.LoS(s_t, avg_amp_gain);
map2d = Map2D();
l4c = Likelihood4Coordinates();
optimiser = gridOptimiser();
ula = ULA(lambda, ELEMENT_NUM, element_spacing);
metric = Metric(); % Initialize the Metric class

% Arrays to store error values for each iteration (for confidence bands)
all_errors = cell(n_param, num_methods);
all_ml_errors = cell(n_param, 1);
for i = 1:n_param
    for j = 1:num_methods
        all_errors{i,j} = zeros(ITERATION, 1);
    end
    all_ml_errors{i} = zeros(ITERATION, 1);
end

rmse_values = zeros(n_param, num_methods);
rmse_values_ml = zeros(n_param, 1);
CRB_values = zeros(n_param, 1);
w = cell(RX_NUM, 1); % Received signal at each Rx vectorised to cell array
aoa_rel_est = zeros(RX_NUM, num_methods);
rays_abs = cell(RX_NUM, num_methods);

%% === Monte Carlo iterations
for itr = 1:ITERATION
    %% --- Location and AoA Refresh for each iteration - ONLY ENABLE FOR RANDOMISED RX AND AOA
    if false % Keep this disabled as we're using fixed positions
        % Generate random RX positions ensuring minimum distance from TX
        pos_rx = zeros(RX_NUM, 2);
        for i = 1:RX_NUM
            valid_position = false;
            while ~valid_position
                % Generate random position
                pos_rx(i,:) = area_size * rand(1, 2);

                % Check if it's far enough from TX
                if sqrt(sum((pos_tx - pos_rx(i,:)).^2)) >= SAFETY_DISTANCE
                    valid_position = true;
                end
            end
        end
        % Generate random true Angle of Arrival from RX to TX
        aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1);
        angle_rx_tx_abs = zeros(RX_NUM, 1);
        for i = 1:RX_NUM
            % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
            angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
        end
        rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
    end
    
    %% === Loop through each SNR value
    for snr_idx=1:n_param
        progressbar('step'); % Update progress bar
        %% === Generate the received signal at each Rx
        for rx_idx=1:RX_NUM
            % --- Generate signal received at Rx
            nPower = avg_E/db2pow(SNR_dB(snr_idx, rx_idx));
            y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
            y_awgn = channel.AWGN(y_ula, nPower);
            % --- append received signal to a centralised array for direct ML estimation
            w{rx_idx} = y_awgn;
        end
        %% === Loop through each method for DoA estimation
        for method_idx = 1:num_methods
            % --- DoA Estimation Algorithm at each RX
            for rx_idx = 1:RX_NUM
                estimator = DoAEstimator(ula, sweeping_angle, aoa_act(rx_idx), DOA_MODE, OPT_GRID_DENSITY);
                aoa_rel_est(rx_idx, method_idx) = estimator.(doa_est_methods(method_idx).name)(y_awgn, doa_est_methods(method_idx).extra_args{:}).aoa_est;
                rays_abs{rx_idx, method_idx} = map2d.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx, method_idx));
            end
            
            % --- Calculate the aoa intersection point and the RMSE for each method
            aoa_intersect = map2d.calDoAIntersect(rays_abs{1, method_idx}, rays_abs{2, method_idx});
            
            % Handle potential intersection failures
            if isempty(aoa_intersect) || ~isfield(aoa_intersect, 'x') || ~isfield(aoa_intersect, 'y')
                error_val = area_size; % Use maximum possible error when rays don't intersect
            else
                error_val = sqrt((pos_tx(1,1)-aoa_intersect.x)^2 + (pos_tx(1,2)-aoa_intersect.y)^2);
            end
            
            % Cap error values using the Metric class
            capped_error = metric.capErrorValues(error_val, area_size/2);
            rmse_values(snr_idx, method_idx) = rmse_values(snr_idx, method_idx) + capped_error;
            all_errors{snr_idx, method_idx}(itr) = capped_error;
        end
        
        %% === direct ML estimation
        objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
        [optCoord, ~] = optimiser.fmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
        
        error_val_ml = sqrt((pos_tx(1,1)-optCoord(1))^2 + (pos_tx(1,2)-optCoord(2))^2);
        capped_error_ml = metric.capErrorValues(error_val_ml, area_size/2);
        rmse_values_ml(snr_idx) = rmse_values_ml(snr_idx) + capped_error_ml;
        all_ml_errors{snr_idx}(itr) = capped_error_ml;
    end
end

% Calculate average RMSE
rmse_values = rmse_values / ITERATION; % RMSE for DoA estimation methods
rmse_values_ml = rmse_values_ml / ITERATION; % RMSE for ML estimation method

%% Calculate confidence bands
lower_bands = zeros(n_param, num_methods);
upper_bands = zeros(n_param, num_methods);

% for method_idx = 1:num_methods
%     for snr_idx = 1:n_param
%         lower_bands(snr_idx, method_idx) = metric.calculateLowerBound(all_errors{snr_idx, method_idx});
%         upper_bands(snr_idx, method_idx) = metric.calculateUpperBound(all_errors{snr_idx, method_idx});
%     end
% end

% % Calculate confidence bands for ML method
% ml_lower_band = zeros(n_param, 1);
% ml_upper_band = zeros(n_param, 1);

% for snr_idx = 1:n_param
%     ml_lower_band(snr_idx) = metric.calculateLowerBound(all_ml_errors{snr_idx});
%     ml_upper_band(snr_idx) = metric.calculateUpperBound(all_ml_errors{snr_idx});
% end

%% === Prepare data for plotting
% Combine RMSE values from DoA methods and ML method for plotting
all_rmse_values = [rmse_values, rmse_values_ml];

% Create display names for all methods
display_names = cell(1, num_methods + 1);
for i = 1:num_methods
    switch DOA_MODE
        case 'sweep'
            extra_str = [' (DoA) ', DOA_MODE, ' with ', num2str(RESOLUTION), ' angle resolution'];
        case 'opt'
            extra_str = [' (DoA) ', DOA_MODE, ' with ', num2str(OPT_GRID_DENSITY), 'x', num2str(OPT_GRID_DENSITY), ' coarse grid'];
        otherwise
            extra_str = '';
    end
    display_names{i} = [strrep(doa_est_methods(i).name, '_', ' '), extra_str];
end
display_names{num_methods + 1} = ['ML coor opt with ', num2str(OPT_GRID_DENSITY), 'x', num2str(OPT_GRID_DENSITY), ' coarse grid'];

% Prepare band data
% all_lower_bands = [lower_bands, ml_lower_band];
% all_upper_bands = [upper_bands, ml_upper_band];
all_lower_bands = [];
all_upper_bands = [];

% Setup show bands flag
show_bands = false(1, num_methods + 1);
if ITERATION > 1
    show_bands = true(1, num_methods + 1);
end

% Use Metric.plots for visualization
metric.plots(mean(SNR_dB, 2), all_rmse_values, 'semilogy', ...
    'DisplayNames', display_names, ...
    'ShowBands', show_bands, ...
    'BandLower', all_lower_bands, ...
    'BandUpper', all_upper_bands, ...
    'Title', ['RMSE of ', num2str(TIME_INST_NUM), 'ts for ', num2str(ITERATION), ' iterations'], ...
    'XLabel', 'Signal to Noise Ratio (SNR) [dB]', ...
    'YLabel', 'Root Mean Square Error (RMSE) [m]', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'LegendLocation', 'northeast');

% Save metrics to file if we have enough iterations for meaningful statistics
% if ITERATION > 1
%     metric_filename = ['metrics_' datestr(now, 'yyyymmdd_HHMMSS') '.mat'];
%     save(['outputs/' metric_filename], 'rmse_values', 'rmse_values_ml', 'lower_bands', 'upper_bands', 'ml_lower_band', 'ml_upper_band', 'SNR_dB');
%     fprintf('Metrics saved to %s\n', metric_filename);
% end
