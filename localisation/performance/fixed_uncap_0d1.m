clear; clc; close all;
%#ok<*UNRCH,*NASGU> % Suppress warnings for unreachable code and unused variables
%% User Inputs and Configurations
ITERATION = 9000;                    % Number of Monte Carlo iterations
RANDOMISE_RX = false;                % Randomise RX positions and AoA
CAP_ERROR = false;                   % Cap error values at the maximum theoretical value
DOA_MODE = 'sweep';                 % DoA estimation mode ('sweep' or 'opt')
DOA_RESOLUTION = 0.1;               % Angle resolution (degrees)
OPT_GRID_DENSITY = 10;              % Define a coarse grid for initial guesses
ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements
NUM_RX_DOA = 2;                     % Number of receivers
SAFETY_DISTANCE = 2;                % Minimum distance between TX and RX (meters)
METRIC_TO_PLOT = 'rmse';            % Options: 'rmse', 'p25', 'p50' (median), 'p75', 'band'
BAND_PERCENTILES = [25, 50, 75];    % Percentiles for error band if METRIC_TO_PLOT is 'band'
SHOW_ERROR_BAND = false;            % Whether to show the 25-75 percentile band
%% Additional RX counts for ML optimization
NUM_RX_ML = 2:8:10;                 % Additional receiver counts for ML optimization
nvar_mlpos = length(NUM_RX_ML);     % Number of variants for ML optimization
%% Initialise classes
channel = ChannelModels();
map2d = Map2D();
metric = Metric();
l4c = Likelihood4Coordinates();
optimiser = Optimisers();
%% Transmitter, receiver positions and angles
area_size = 100;
pos_tx = [50, 50];
%% SNR values to test
SNR_dB = repmat((-10:2:20)', 1, max(NUM_RX_DOA, max(NUM_RX_ML)));    % SNR in dB
nvar_snr = length(SNR_dB);                   % Number of positions to test
%% Signal and channel configurations
c = 299792458;                              % Speed of light (m/s)
fc = 2.4e9;                                 % Base carrier frequency (Hz) (known)
lambda = c / fc;                            % Wavelength (m)
avg_amp_gain = 1;                           % Average gain of the channel
L_d0=100;                                   % Reference Power (dB) - for gain calculation
d0=100;                                     % Reference distance (m) - for gain calculation
alpha=4;                                    % Path loss exponent - for gain calculation
P_t = 1;                                    % W - Transmit signal power (known)                      
Fs = 2 * fc;                                % Sample frequency, enough for the signal
T = TIME_INST_NUM/Fs;                       % Period of transmission
t = 0:1/Fs:(T-1/Fs);                        % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;             % Element spacing (ULA)
sweeping_angle = -90:DOA_RESOLUTION:90;         % Angle range for finding the AoA

% Generate nuisance transmitted signal with random phase
% Only the frequency and power are known
rng('shuffle');                             % Ensure randomness on each run
random_phases = 2*pi*rand(1, TIME_INST_NUM); % Random phase for each time instance
random_amp_variations = 0.2*randn(1, TIME_INST_NUM) + 1; % Minor amplitude variations

% Construct the nuisance signal
s_t = zeros(1, length(t));
for i = 1:TIME_INST_NUM
    % Get sample indices for this time instance
    start_idx = round((i-1)*length(t)/TIME_INST_NUM) + 1;
    end_idx = round(i*length(t)/TIME_INST_NUM);
    
    % Create signal with known frequency and power but random phase
    s_t(start_idx:end_idx) = sqrt(P_t) * random_amp_variations(i) * ...
        exp(1j * (2*pi*fc*t(start_idx:end_idx) + random_phases(i)));
end

% Calculate average energy of the signal
e_avg = FIXED_TRANS_ENERGY * 1 + ~FIXED_TRANS_ENERGY * mean(abs(s_t).^2) * T;
%% === Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'BF'}, ... % estimator methods
    'extra_args', {{}} ...        % extra args required for specific type of estimator
);
nvar_doa = numel(doa_est_methods);               % Automatically get number of methods from struct array
num_legend = nvar_doa + nvar_mlpos;         % Number of methods plus ML method
fprintf('Running Monte Carlo simulation with %d iterations...\n', ITERATION);
progressbar('reset', ITERATION*nvar_snr+ ITERATION*nvar_mlpos*nvar_snr);            % Reset progress bar
%% Initialise arrays
ula = ULA(lambda, ELEMENT_NUM, element_spacing);    % Create Uniform Linear Array object

% Initialize storage for all individual errors
all_errors = cell(nvar_snr, num_legend);             % Pre-allocate for variable-sized collections
for i=1:nvar_snr
    for j=1:num_legend
        all_errors{i,j} = zeros(ITERATION, 1);
    end
end
aoa_rel_est = zeros(NUM_RX_DOA, nvar_doa);           % Pre-allocate for Relative AoA estimation
rays_abs = cell(NUM_RX_DOA, nvar_doa);               % Pre-allocate for absolute rays
%% === Monte Carlo iterations
for itr = 1:ITERATION
    %% --- DoA Estimation
    [pos_rx, aoa_act, rot_abs] = map2d.genPos(area_size, pos_tx, NUM_RX_DOA, RANDOMISE_RX, SAFETY_DISTANCE, ABS_ANGLE_LIM, DOA_RESOLUTION);
    [~, y_centralised] = channel.generateReceivedSignal(s_t, pos_tx, pos_rx, aoa_act, e_avg, SNR_dB, L_d0, d0, alpha, ELEMENT_NUM, element_spacing, lambda);
    for idx_snr=1:nvar_snr % Loop through each SNR value
        progressbar('step'); % Update progress bar
        for idx_doa_method = 1:nvar_doa
            % --- Estimate the AoA and the ray to that AoA for each receiver
            for idx_rx = 1:NUM_RX_DOA
                estimator = DoAEstimator(ula, sweeping_angle, aoa_act(idx_rx), DOA_MODE, OPT_GRID_DENSITY);
                aoa_rel_est(idx_rx, idx_doa_method) = estimator.(doa_est_methods(idx_doa_method).name)(y_centralised{idx_snr, idx_rx}, doa_est_methods(idx_doa_method).extra_args{:}).aoa_est;
                rays_abs{idx_rx, idx_doa_method} = map2d.calAbsRays(pos_rx(idx_rx,:), pos_tx, rot_abs(idx_rx), aoa_rel_est(idx_rx, idx_doa_method));
            end
            % --- Calculate the aoa intersection point and the error distance
            aoa_intersect = map2d.calDoAIntersect(rays_abs{1, idx_doa_method}, rays_abs{2, idx_doa_method});
            all_errors{idx_snr, idx_doa_method}(itr) = sqrt((pos_tx(1,1)-aoa_intersect.x)^2 + (pos_tx(1,2)-aoa_intersect.y)^2);
        end
    end
    %% --- ML optimization with additional receivers
    for ml_idx = 1:nvar_mlpos
        % --- Generate receivers and the received signal
        [pos_rx_ml, aoa_act_ml, rot_abs_ml] = map2d.genPos(area_size, pos_tx, NUM_RX_ML(ml_idx), RANDOMISE_RX, SAFETY_DISTANCE, ABS_ANGLE_LIM, DOA_RESOLUTION);
        [nPower, y_centralised_ml] = channel.generateReceivedSignal(s_t, pos_tx, pos_rx_ml, aoa_act_ml, e_avg, SNR_dB, L_d0, d0, alpha, ELEMENT_NUM, element_spacing, lambda);
        % Loop through each SNR value
        for idx_snr=1:nvar_snr
            progressbar('step'); % Update progress bar
            % % --- Direct ML estimation
            objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx_ml, rot_abs_ml, y_centralised_ml(idx_snr, :)', ELEMENT_NUM, nPower);
            % [optCoord, ~] = optimiser.gridFmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
            % all_errors{idx_snr, nvar_doa+ml_idx}(itr) = sqrt((pos_tx(1,1)-optCoord(1))^2 + (pos_tx(1,2)-optCoord(2))^2);
            % --- Estimate the AoA and the ray to that AoA for each receiver
            for idx_rx = 1:NUM_RX_DOA
                estimator = DoAEstimator(ula, sweeping_angle, aoa_act_ml(idx_rx), DOA_MODE, OPT_GRID_DENSITY);
                aoa_rel_est(idx_rx, 1) = estimator.(doa_est_methods(1).name)(y_centralised_ml{idx_snr, idx_rx}, doa_est_methods(1).extra_args{:}).aoa_est;
                rays_abs{idx_rx, 1} = map2d.calAbsRays(pos_rx_ml(idx_rx,:), pos_tx, rot_abs_ml(idx_rx), aoa_rel_est(idx_rx, 1));
            end
            % --- Calculate the aoa intersection point and the error distance
            aoa_intersect = map2d.calDoAIntersect(rays_abs{1, 1}, rays_abs{2, 1});
            [optCoord, ~] = optimiser.findMinWithInitialPoint(objective_to_maximize, {}, [0, 0], [area_size, area_size], [aoa_intersect.x, aoa_intersect.y]);
            all_errors{idx_snr, nvar_doa+ml_idx}(itr) = sqrt((pos_tx(1,1)-optCoord(1))^2 + (pos_tx(1,2)-optCoord(2))^2);
        end
    end
end
% Cap errors at the maximum theoretical value
if CAP_ERROR
    max_possible_error = sqrt(2) * area_size;
    all_errors = metric.capErrorValues(all_errors, max_possible_error);
end

percentiles = struct( ...
    'lower', zeros(nvar_snr, num_legend), ...
    'upper', zeros(nvar_snr, num_legend), ...
    'val', zeros(nvar_snr, num_legend));
% Select which metric to calculate and plot
switch METRIC_TO_PLOT
    case 'rmse'
        plot_data = metric.cal_RMSE(all_errors);
        metric_label = 'RMSE';
    case 'p25'
        plot_data = metric.cal_Percentiles(all_errors, 25).val;
        metric_label = '25th Percentile';
    case 'p75'
        plot_data = metric.cal_Percentiles(all_errors, 75).val;
        metric_label = '75th Percentile';
    case 'band' % Calculate the percentiles for the error band
        percentiles = metric.cal_Percentiles(all_errors, BAND_PERCENTILES);
        plot_data = percentiles.val;
        metric_label = 'Median with Error Band';
    otherwise % Default to median (50th percentile)
        METRIC_TO_PLOT = 'p50';
        plot_data = metric.cal_Percentiles(all_errors).val;
        metric_label = 'Median';
end

%% === Prepare data and plotting
% Create display names for all methods
legend_name = cell(1, num_legend);
for i = 1:nvar_doa
    switch DOA_MODE
        case 'sweep'
            modeString = [DOA_MODE, ' (', num2str(DOA_RESOLUTION), '\circ res)'];
        case 'opt'
            modeString = [DOA_MODE, ' (', num2str(OPT_GRID_DENSITY), ' grid)'];
        otherwise
            modeString = [DOA_MODE];
    end
    legend_name{i} = [strrep(doa_est_methods(i).name, '_', ' '), ' DoA for ', num2str(NUM_RX_DOA), ' RXs by ', modeString];
end
for ml_idx = 1:nvar_mlpos
    legend_name{nvar_doa+ml_idx} = ['MLpos of ' num2str(NUM_RX_ML(ml_idx)) ' RXs from 2 DoA initial'];
end
rx_type = {'fixed', 'randomised'};
cap_error = {'uncapped', 'capped'};
annotStrings = {
    ['RX Type: ', rx_type{1 + RANDOMISE_RX}], ...
    ['ULA elements: ', num2str(ELEMENT_NUM)], ...
    ['Time instances: ', num2str(TIME_INST_NUM)], ...
    ['Error Metric: ', metric_label, ' (', cap_error{1+CAP_ERROR},')']
};

% Plot the error metric
metric.plots(mean(SNR_dB, 2), plot_data, 'semilogy', ...
    'DisplayNames', legend_name, ...
    'ShowBands', strcmp(METRIC_TO_PLOT, 'band') * ones(1, num_legend), ...
    'BandLower', percentiles.lower, ...
    'BandUpper', percentiles.upper, ...
    'Title', ['Error Metric by estimation method (', num2str(ITERATION), ' iterations)'], ...
    'YLabel', [metric_label, ' Error [m]'], ...
    'ShowAnnotation', true, ...
    'AnnotationStrings', annotStrings);

%% Save metrics to file if we have enough iterations for meaningful statistics
if ITERATION > 0
OUTPUT_PATH = 'data';
    % Create the outputs directory if it doesn't exist
    if ~exist(OUTPUT_PATH, 'dir')
        mkdir(OUTPUT_PATH);
    end
    
    % Generate filename with datetime in standard format
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    metric_filename = [timestamp '_' mfilename '.mat'];
    
    % Save the data including annotation strings and display names
    save( ...
        fullfile(OUTPUT_PATH, metric_filename), ...
        'SNR_dB', 'plot_data', 'all_errors', 'ITERATION', ...
        'legend_name', 'METRIC_TO_PLOT', 'num_legend', 'percentiles', ...
        'metric_label', 'annotStrings');
    fprintf('Metrics saved to %s\n', metric_filename);
end
