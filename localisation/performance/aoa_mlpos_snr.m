clear; clc; close all;
%#ok<*UNRCH,*NASGU> % Suppress warnings for unreachable code and unused variables

%% User Inputs and Configurations
RUN_MODE = 'test'; % Options: 'test' or 'full'

% Set ITERATION based on run mode
if strcmp(RUN_MODE, 'test')
    ITERATION = 1;  % For quick testing
    SAVE_METRICS = false;
    SHOW_PLOTS = true;
else  % 'full' mode
    ITERATION = 30000;  % For full simulation
    SAVE_METRICS = true;
    SHOW_PLOTS = false;
end

RANDOMISE_RX = false;               % Randomise RX positions and AoA
CAP_ERROR = false;                   % Cap error values at the maximum theoretical value
INCLUDE_CAPPED = true;             % Include capped values in the output errors, only valid if CAP_ERROR is true
COMPARE_EPDF_IN_SUBPLOT = true;     % Compare empirical PDFs in subplots
DOA_MODE = 'sweep';                 % DoA estimation mode ('sweep' or 'opt')
DOA_RESOLUTION = 1;               % Angle resolution (degrees)
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
NUM_RX_ML = 3:7:10;                 % Additional receiver counts for ML optimization
nvar_mlpos = length(NUM_RX_ML);     % Number of variants for ML optimization

%% Initialize classes
channel = ChannelModels();
metric = Metric();
l4c = Likelihood4Coordinates();
optimiser = Optimisers();
algo = Algorithms(l4c, optimiser);
map2d = Map2D([10,10], [90, 90], max(NUM_RX_ML));

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
if ~CAP_ERROR
    INCLUDE_CAPPED = true;                 % Disable capped values if we're not capping errors
end

% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;             % Element spacing (ULA)
sweeping_angle = -90:DOA_RESOLUTION:90;     % Angle range for finding the AoA

% Generate nuisance transmitted signal with random phase
[s_t, e_avg] = channel.generateNuisanceSignal(fc, P_t, T, t, TIME_INST_NUM, FIXED_TRANS_ENERGY);

%% === Define the methods to test for performance
doa_est_methods = struct(...
    'name', {'BF'}, ... % estimator methods
    'extra_args', {{}} ...        % extra args required for specific type of estimator
);
nvar_doa = numel(doa_est_methods);               % Automatically get number of methods from struct array
num_legend = nvar_doa + nvar_mlpos*2;         % Number of methods plus ML method

fprintf('Running Monte Carlo simulation with %d iterations...\n', ITERATION);
progressbar('reset', ITERATION*nvar_snr+ ITERATION*nvar_mlpos*2*nvar_snr);            % Reset progress bar

%% Initialise arrays
ula = ULA(lambda, ELEMENT_NUM, element_spacing);    % Create Uniform Linear Array object
estimator = DoAEstimator(ula, sweeping_angle, 0, DOA_MODE, OPT_GRID_DENSITY);

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
    % --- Generate receivers and the received signal for the maximum number of receivers
    [pos_rx, aoa_act, rot_abs] = map2d.genRXPos(area_size, pos_tx, max(NUM_RX_ML), RANDOMISE_RX, SAFETY_DISTANCE, ABS_ANGLE_LIM, DOA_RESOLUTION);
    [nPower, y_centralised] = channel.generateReceivedSignal(s_t, pos_tx, pos_rx, aoa_act, e_avg, SNR_dB, L_d0, d0, alpha, ELEMENT_NUM, element_spacing, lambda);
    doa_estimator = @(sig) estimator.(doa_est_methods(1).name)(sig, doa_est_methods(1).extra_args{:});
    
    %% --- DoA estimation only
    for idx_snr=1:nvar_snr
        y_received = y_centralised(idx_snr, :);
        progressbar('step'); % Update progress bar
        [~, all_errors{idx_snr, 1}(itr)] = algo.DoAintersect(...
            pos_rx, rot_abs, y_received, ...
            doa_estimator, pos_tx ...
        );
    end
    
    %% --- ML optimization with additional receivers
    for ml_idx = 1:nvar_mlpos
        % Loop through each SNR value
        for idx_snr=1:nvar_snr
            pos_rx_active = pos_rx(1:NUM_RX_ML(ml_idx),:);
            y_received_active = y_centralised(idx_snr, 1:NUM_RX_ML(ml_idx),:);
            progressbar('step'); % Update progress bar
            [~, ~, all_errors{idx_snr, nvar_doa+ml_idx}(itr)] = algo.MLOpt4mDoAtriage(...
                pos_rx_active, rot_abs, y_received_active, ...
                ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                doa_estimator, pos_tx...
            );
            progressbar('step'); % Update progress bar
            [~, ~, all_errors{idx_snr, nvar_doa+nvar_mlpos+ml_idx}(itr)] = algo.MLOpt4mCentroid(...
                pos_rx_active, rot_abs, y_received_active, ...
                ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                doa_estimator, pos_tx...
            );
            % progressbar('step'); % Update progress bar
            % [~, ~, all_errors{idx_snr, nvar_doa+2*nvar_mlpos+ml_idx}(itr)] = algo.MLOptwGrid(...
            %     pos_rx_active, rot_abs, y_received_active, ...
            %     ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
            %     OPT_GRID_DENSITY, pos_tx...
            % );
        end
    end
end

%% Process results - either save or plot
if SAVE_METRICS
    % Save metrics to file
    OUTPUT_PATH = 'data';
    % Create the outputs directory if it doesn't exist
    if ~exist(OUTPUT_PATH, 'dir')
        mkdir(OUTPUT_PATH);
    end
    
    % Generate filename with datetime in standard format
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    metric_filename = [timestamp '_' mfilename '.mat'];
    
    % Save the data including all parameters needed for plotting
    save(fullfile(OUTPUT_PATH, metric_filename));
    fprintf('Metrics saved to data/%s\n', metric_filename);
end

if SHOW_PLOTS
    % Call the shared plotting function for immediate visualization
    plot_aoa_mlpos(all_errors, SNR_dB, metric, doa_est_methods, ...
        NUM_RX_DOA, NUM_RX_ML, DOA_MODE, DOA_RESOLUTION, OPT_GRID_DENSITY, ...
        ELEMENT_NUM, TIME_INST_NUM, RANDOMISE_RX, CAP_ERROR, INCLUDE_CAPPED, ...
        COMPARE_EPDF_IN_SUBPLOT, ITERATION, METRIC_TO_PLOT, BAND_PERCENTILES);
end
