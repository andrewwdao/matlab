%#ok<NASGU>
%#ok<UNRCH>
clear; clc; close all;
%% User Inputs and Configurations
ITERATION = 1500; % Number of Monte Carlo iterations
OPT_GRID_DENSITY = 10; % Define a coarse grid for initial guesses
ABS_ANGLE_LIM = 60;                  % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 0.1;                   % Angle resolution (degrees)
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements
DOA_MODE = 'sweep';                 % DoA estimation mode ('sweep' or 'opt')
SAFETY_DISTANCE = 2;                % Minimum distance between TX and RX (meters)
RANDOMISE_RX = true;               % Randomise RX positions and AoA
RX_NUM = 2;                        % Number of receivers
SHOW_ERROR_BAND = false;            % Whether to show the 25-75 percentile band
% Add parameter to select which metric to plot
METRIC_TO_PLOT = 'p50';             % Options: 'rmse', 'p25', 'p50' (median), 'p75', 'band'
BAND_PERCENTILES = [25, 50, 75];        % Percentiles for error band if METRIC_TO_PLOT is 'band'

%% Initialise classes
channel = ChannelModels();
map2d = Map2D();
metric = Metric();
l4c = Likelihood4Coordinates();
optimiser = gridOptimiser();
%% Transmitter, receiver positions and angles
area_size = 100;
pos_tx = [50, 50];
if ~RANDOMISE_RX % Fixed rx and aoa
    % pos_rx = [21, 51; 21, 60]; % 1
    % pos_rx = [21, 51; 30, 70]; % 2
    % pos_rx = [21, 51; 40, 70]; % 3
    % pos_rx = [21, 51; 50, 70]; % 4
    % pos_rx = [21, 51; 60, 70]; % 5
    % pos_rx = [21, 51; 70, 60]; % 6
    % pos_rx = [21, 51; 70, 50]; % 7
    % pos_rx = [21, 51; 70, 40]; % 8
    % pos_rx = [21, 51; 60, 30]; % 9
    % pos_rx = [21, 51; 50, 30]; % 10
    % pos_rx = [21, 51; 40, 30]; % 11
    % pos_rx = [21, 51; 30, 30]; % 12
    pos_rx = [21, 51; 20, 40]; % 2 Rx at fixed positions
    aoa_act = [0; 0];            % True AoA from Rx to Tx
    RX_NUM = size(pos_rx, 1);    % Update RX_NUM based on number of receivers
    rot_abs = map2d.calAbsAngle(pos_tx, pos_rx, aoa_act);
end

% --- Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
% angle_rx_tx_abs = zeros(RX_NUM, 1);
% for i = 1:RX_NUM
%     angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
% end
% rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees

%% SNR values to test
SNR_dB = repmat((-10:2:20)', 1, RX_NUM);       % SNR in dB
n_param = length(SNR_dB); % Number of positions to test
%% Signal and channel configurations
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(RX_NUM, 1);  % W - Transmit signal power
sub_carrier = (1:RX_NUM)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA

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
num_legend = num_methods + 1; % Number of methods plus ML method
fprintf('Running Monte Carlo simulation with %d iterations...\n', ITERATION);
progressbar('reset', ITERATION*n_param); % Reset progress bar
%% Initialise arrays
y_los = channel.LoS(s_t, avg_amp_gain);
ula = ULA(lambda, ELEMENT_NUM, element_spacing);

% Arrays to store error values for each iteration
all_errors = cell(n_param, num_methods);
all_ml_errors = cell(n_param, 1);
for i = 1:n_param
    for j = 1:num_methods
        all_errors{i,j} = zeros(ITERATION, 1);
    end
    all_ml_errors{i} = zeros(ITERATION, 1);
end

centralised_received_signal = cell(RX_NUM, 1); % Received signal at each Rx vectorised to cell array
% rmse_values = zeros(n_param, num_methods);
% rmse_values_ml = zeros(n_param, 1);
aoa_rel_est = zeros(RX_NUM, num_methods);
rays_abs = cell(RX_NUM, num_methods);
%% === Monte Carlo iterations
for itr = 1:ITERATION
    %% --- Location and AoA Refresh for each iteration - ONLY ENABLE FOR RANDOMISED RX AND AOA
    if RANDOMISE_RX
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
        rot_abs = map2d.calAbsAngle(pos_tx, pos_rx, aoa_act);
        % angle_rx_tx_abs = zeros(RX_NUM, 1);
        % for i = 1:RX_NUM
        %     % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
        %     angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
        % end
        % rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
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
            centralised_received_signal{rx_idx} = y_awgn;
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
            
            % % Handle potential intersection failures
            % if isempty(aoa_intersect) || ~isfield(aoa_intersect, 'x') || ~isfield(aoa_intersect, 'y')
            %     error_val = area_size; % Use maximum possible error when rays don't intersect
            % else
            %     error_val = sqrt((pos_tx(1,1)-aoa_intersect.x)^2 + (pos_tx(1,2)-aoa_intersect.y)^2);
            % end
            all_errors{snr_idx, method_idx}(itr) = sqrt((pos_tx(1,1)-aoa_intersect.x)^2 + (pos_tx(1,2)-aoa_intersect.y)^2);
            % Cap error values using the Metric class
            % capped_error = metric.capErrorValues(error_val, area_size/2);
            % rmse_values(snr_idx, method_idx) = rmse_values(snr_idx, method_idx) + capped_error;
            % all_errors{snr_idx, method_idx}(itr) = capped_error;
        end
        
        %% === direct ML estimation
        objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, centralised_received_signal, ELEMENT_NUM, nPower);
        [optCoord, ~] = optimiser.fmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
        all_ml_errors{snr_idx}(itr) = sqrt((pos_tx(1,1)-optCoord(1))^2 + (pos_tx(1,2)-optCoord(2))^2);
        % error_val_ml = sqrt((pos_tx(1,1)-optCoord(1))^2 + (pos_tx(1,2)-optCoord(2))^2);
        % capped_error_ml = metric.capErrorValues(error_val_ml, area_size/2);
        % rmse_values_ml(snr_idx) = rmse_values_ml(snr_idx) + capped_error_ml;
        % all_ml_errors{snr_idx}(itr) = capped_error_ml;
    end
end
% Cap errors at the maximum theoretical value
max_possible_error = sqrt(2) * area_size;
all_errors = metric.capErrorValues(all_errors, max_possible_error);
all_ml_errors = metric.capErrorValues(all_ml_errors, max_possible_error);
% Calculate average RMSE
% rmse_values = rmse_values / ITERATION; % RMSE for DoA estimation methods
% rmse_values_ml = rmse_values_ml / ITERATION; % RMSE for ML estimation method

percentiles = struct( ...
    'lower', zeros(n_param, num_legend), ...
    'upper', zeros(n_param, num_legend), ...
    'val', zeros(n_param, num_legend));
% Select which metric to calculate and plot
switch METRIC_TO_PLOT
    case 'rmse'
        plot_data = [
            metric.cal_RMSE(all_errors), ...
            metric.cal_RMSE(all_ml_errors)
        ];
        metric_label = 'RMSE';
    case 'p25'
        plot_data = [
            metric.cal_Percentiles(all_errors, 25).val, ...
            metric.cal_Percentiles(all_ml_errors, 25).val
        ];
        metric_label = '25th Percentile';
    case 'p75'
        plot_data = [
            metric.cal_Percentiles(all_errors, 75).val, ...
            metric.cal_Percentiles(all_ml_errors, 75).val
        ];
        metric_label = '75th Percentile';
    case 'band'
        % Calculate the 25th, 50th, and 75th percentiles
        percentiles = [
            metric.cal_Percentiles(all_errors, BAND_PERCENTILES), ...
            metric.cal_Percentiles(all_ml_errors, BAND_PERCENTILES)
        ];
        % Plot the median with error band
        plot_data = percentiles.value;
        metric_label = 'Median with Error Band';
    otherwise
        % Default to median (50th percentile)
        METRIC_TO_PLOT = 'p50';
        plot_data = [
            metric.cal_Percentiles(all_errors).val, ...
            metric.cal_Percentiles(all_ml_errors).val
        ];
        metric_label = 'Median';
end

% %% Calculate confidence bands
% lower_bands = zeros(n_param, num_methods);
% upper_bands = zeros(n_param, num_methods);

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
% all_rmse_values = [rmse_values, metric.cal_RMSE(all_errors), rmse_values_ml, metric.cal_RMSE(all_ml_errors)];

% Create display names for all methods
display_names = cell(1, num_methods + 3);
for i = 1:num_methods
    switch DOA_MODE
        case 'sweep'
            % extra_str = [' (DoA) ', DOA_MODE, ' with ', num2str(RESOLUTION), ' angle resolution'];
            modeString = ['Mode: ', DOA_MODE, ' (', num2str(RESOLUTION), '\circ res)'];
        case 'opt'
            % extra_str = [' (DoA) ', DOA_MODE, ' with ', num2str(OPT_GRID_DENSITY), 'x', num2str(OPT_GRID_DENSITY), ' coarse grid'];
            modeString = ['Mode: ', DOA_MODE, ' (', num2str(OPT_GRID_DENSITY), ' grid)'];
        otherwise
            % extra_str = '';
            modeString = ['Mode: ', DOA_MODE];
    end
    display_names{i} = [strrep(doa_est_methods(i).name, '_', ' '), ' (DoA) '];
end
display_names{num_legend} = ['ML coor opt with ', num2str(OPT_GRID_DENSITY), 'x', num2str(OPT_GRID_DENSITY), ' grid'];
annotStrings = {
    ['RX number: ', num2str(RX_NUM)], ...    
    ['ULA elements: ', num2str(ELEMENT_NUM)], ...
    ['Time instances: ', num2str(TIME_INST_NUM)], ...
    modeString, ...
    ['Error Metric: ', metric_label]
};
% % Prepare band data
% % all_lower_bands = [lower_bands, ml_lower_band];
% % all_upper_bands = [upper_bands, ml_upper_band];
% all_lower_bands = [];
% all_upper_bands = [];

% % Setup show bands flag
% show_bands = false(1, num_methods + 3);
% % if ITERATION > 1
% %     show_bands = true(1, num_methods + 3);
% % end

% Use Metric.plots for visualization
metric.plots(mean(SNR_dB, 2), plot_data, 'semilogy', ...
    'DisplayNames', display_names, ...
    'ShowBands', strcmp(METRIC_TO_PLOT, 'band') * ones(1, num_legend), ...
    'BandLower', percentiles.lower, ...
    'BandUpper', percentiles.upper, ...
    'Title', ['Error Metric (Cap) by estimation method (', num2str(ITERATION), ' iterations)'], ...
    'YLabel', [metric_label, ' Error [m]'], ...
    'ShowAnnotation', true, ...
    'AnnotationStrings', annotStrings);

% Save metrics to file if we have enough iterations for meaningful statistics
% if ITERATION > 1
%     metric_filename = ['metrics_' datestr(now, 'yyyymmdd_HHMMSS') '.mat'];
%     save(['outputs/' metric_filename], 'rmse_values', 'rmse_values_ml', 'lower_bands', 'upper_bands', 'ml_lower_band', 'ml_upper_band', 'SNR_dB');
%     fprintf('Metrics saved to %s\n', metric_filename);
% end
