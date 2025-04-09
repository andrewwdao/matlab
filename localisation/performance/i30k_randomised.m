%#ok<*UNRCH,*NASGU> % Suppress warnings for unreachable code and unused variables
clear; clc; close all;

%% User Inputs and Configurations
RUN_MODE = 'save';                  % Options: 'test', 'plot' or 'save'
METRIC_TO_PLOT = 'rmse';            % Options: 'rmse', 'p25', 'p50' (median), 'p75', 'band'
BAND_PERCENTILES = [25, 50, 75];    % Percentiles for error band if METRIC_TO_PLOT is 'band'
SHOW_ERROR_BAND = false;            % Whether to show the 25-75 percentile band
CAP_ERROR = false;                  % Cap error values at the maximum theoretical value
INCLUDE_CAPPED = false;              % Include capped values in the output errors, only valid if CAP_ERROR is true
COMPARE_EPDF_IN_SUBPLOT = false;     % Compare empirical PDFs in subplots

FLAG_PLOT = strcmp(RUN_MODE, 'test') || strcmp(RUN_MODE, 'plot'); % Flag for plotting
if strcmp(RUN_MODE, 'plot')
    % Save current UI settings before loading
    original_settings = struct(...
        'METRIC_TO_PLOT', METRIC_TO_PLOT, ...
        'BAND_PERCENTILES', BAND_PERCENTILES, ...
        'SHOW_ERROR_BAND', SHOW_ERROR_BAND, ...
        'CAP_ERROR', CAP_ERROR, ...
        'INCLUDE_CAPPED', INCLUDE_CAPPED, ...
        'COMPARE_EPDF_IN_SUBPLOT', COMPARE_EPDF_IN_SUBPLOT, ...
        'FLAG_PLOT', FLAG_PLOT);

    % Select which data file to load
    [filename, pathname] = uigetfile('data/*.mat', 'Select a saved simulation result');
    if isequal(filename, 0)
        error('No file selected');
    end
    
    % Load the data file
    load(fullfile(pathname, filename));
    
    % Restore original UI settings
    METRIC_TO_PLOT = original_settings.METRIC_TO_PLOT;
    BAND_PERCENTILES = original_settings.BAND_PERCENTILES;
    SHOW_ERROR_BAND = original_settings.SHOW_ERROR_BAND;
    CAP_ERROR = original_settings.CAP_ERROR;
    INCLUDE_CAPPED = original_settings.INCLUDE_CAPPED;
    COMPARE_EPDF_IN_SUBPLOT = original_settings.COMPARE_EPDF_IN_SUBPLOT;
    FLAG_PLOT = original_settings.FLAG_PLOT;
else
    % Set ITERATION based on run mode
    if strcmp(RUN_MODE, 'test')
        ITERATION = 1;  % For quick testing
        SAVE_METRICS = false;
    else  % 'full' mode
        ITERATION = 30000;  % For full simulation
        SAVE_METRICS = true;
    end

    TX_RANDOMISED = false;              % Randomise TX positions
    RX_RANDOMISED = true;              % Randomise RX positions and AoA
    TX_NUM = 1;                         % Number of transmitters
    RX_NUM = 3:7:24;                 % Additional receiver counts for ML optimization
    nvar_mlpos = length(RX_NUM);     % Number of variants for ML optimization
    ELEMENT_NUM = 4;                    % Number of ULA elements
    DOA_MODE = 'sweep';                 % DoA estimation mode ('sweep' or 'opt')
    DOA_RESOLUTION = 1;                 % Angle resolution (degrees)
    ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
    TIME_INST_NUM = 1;                  % Number of time instances
    FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
    OPT_GRID_DENSITY = 5;              % Define a coarse grid for initial guesses
    SAFETY_DISTANCE = 2;                % Minimum distance between TX and RX (meters)
    area_size = 100;                    % Area size for RX positions
    

    %% Signal and channel configurations
    SNR_dB = repmat((-10:2:20)', 1, max(2, max(RX_NUM)));    % SNR in dB
    nvar_snr = size(SNR_dB, 1);                   % Number of positions to test
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
    sweeping_angle = -90:DOA_RESOLUTION:90;     % Angle range for finding the AoA
    if ~CAP_ERROR
        INCLUDE_CAPPED = true;                 % Disable capped values if we're not capping errors
    end

    %% Initialize classes
    channel = ChannelModels();
    l4c = Likelihood4Coordinates();
    optimiser = Optimisers();
    algo = Algorithms(l4c, optimiser);
    map2d = Map2D([10,10], [90, 90], max(RX_NUM));
    metric = Metric();

    %% The methods to test for performance
    doa_est_methods = struct(...
        'name', {'BF'}, ... % estimator methods
        'extra_args', {{}} ...        % extra args required for specific type of estimator
    );
    nvar_doa = numel(doa_est_methods);               % Automatically get number of methods from struct array
    num_legend = nvar_doa + nvar_mlpos*2;         % Number of methods plus ML method

    % Initialize storage for all individual errors
    all_errors = cell(nvar_snr, num_legend);             % Pre-allocate for variable-sized collections
    for i=1:nvar_snr
        for j=1:num_legend
            all_errors{i,j} = zeros(ITERATION, 1);
        end
    end

    %% Transmitter, receiver positions and angles
    % Transmitters
    pos_tx = map2d.genTXPos(area_size, TX_NUM, TX_RANDOMISED);
    [s_t, e_avg] = channel.generateNuisanceSignal(fc, P_t, T, t, TIME_INST_NUM, FIXED_TRANS_ENERGY); % Generate nuisance transmitted signal with random phase
    % Initialise DoA estimator
    ula = ULA(lambda, ELEMENT_NUM, element_spacing);    % Create Uniform Linear Array object
    estimator = DoAEstimator(ula, sweeping_angle, 0, DOA_MODE, OPT_GRID_DENSITY);
    doa_estimator = @(sig) estimator.(doa_est_methods(1).name)(sig, doa_est_methods(1).extra_args{:});

    fprintf('Running Monte Carlo simulation with %d iterations...\n', ITERATION);
    progressbar('reset', ITERATION*nvar_snr+ ITERATION*nvar_mlpos*2*nvar_snr); % Reset progress bar
    %% === Monte Carlo iterations
    for itr = 1:ITERATION
        % --- Generate receivers and the received signal for the maximum number of receivers
        [pos_rx, aoa_act, rot_abs] = map2d.genRXPos(area_size, pos_tx, max(RX_NUM), RX_RANDOMISED, SAFETY_DISTANCE, ABS_ANGLE_LIM, DOA_RESOLUTION);
        [nPower, y_centralised] = channel.generateReceivedSignal(s_t, pos_tx, pos_rx, aoa_act, e_avg, SNR_dB, L_d0, d0, alpha, ELEMENT_NUM, element_spacing, lambda);

        %% --- DoA estimation only
        for idx_snr=1:nvar_snr
            y_received = y_centralised(idx_snr, :);
            progressbar('step'); % Update progress bar
            % Time the DoAtriage algorithm
            progressbar('starttimer', 'DoA_triage');
            [~, all_errors{idx_snr, 1}(itr)] = algo.DoAtriage(...
                pos_rx, rot_abs, y_received, ...
                doa_estimator, pos_tx ...
            );
            progressbar('stoptimer', 'DoA_triage');
        end
        
        %% --- ML optimization with additional receivers
        for ml_idx = 1:nvar_mlpos
            % Create algorithm names with RX count
            triage_timer_name = sprintf('triage_initial/ML%dRXs', RX_NUM(ml_idx));
            centroid_timer_name = sprintf('centroid_initial/ML%dRXs', RX_NUM(ml_idx));
            % Loop through each SNR value
            for idx_snr=1:nvar_snr
                pos_rx_active = pos_rx(1:RX_NUM(ml_idx),:);
                y_received_active = y_centralised(idx_snr, 1:RX_NUM(ml_idx),:);
                progressbar('step'); % Update progress bar
                progressbar('starttimer', triage_timer_name); % Time the algorithm
                [~, ~, all_errors{idx_snr, nvar_doa+ml_idx}(itr)] = algo.MLOpt4mDoAtriage(...
                    pos_rx_active, rot_abs, y_received_active, ...
                    ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                    doa_estimator, pos_tx...
                );
                progressbar('stoptimer', triage_timer_name);
                progressbar('step'); % Update progress bar
                progressbar('starttimer', centroid_timer_name); % Time the algorithm
                [~, ~, ~, all_errors{idx_snr, nvar_doa+nvar_mlpos+ml_idx}(itr)] = algo.MLOpt4mCentroid(...
                    pos_rx_active, rot_abs, y_received_active, ...
                    ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                    doa_estimator, pos_tx...
                );
                progressbar('stoptimer', centroid_timer_name);
                % progressbar('step'); % Update progress bar
                % [~, ~, all_errors{idx_snr, nvar_doa+2*nvar_mlpos+ml_idx}(itr)] = algo.MLOptwGrid(...
                %     pos_rx_active, rot_abs, y_received_active, ...
                %     ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                %     OPT_GRID_DENSITY, pos_tx...
                % );
            end
        end
    end
% After all iterations, report the execution times
timer_stats = progressbar('reporttimers');
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
end

if FLAG_PLOT
    fprintf('Plotting results...\n');
    nvar_snr = size(SNR_dB, 1);
    nvar_doa = numel(doa_est_methods);
    nvar_mlpos = length(RX_NUM);
    num_legend = nvar_doa + nvar_mlpos*2;
    
    %% Cap errors at the maximum theoretical value if needed
    area_size = 100; % This should be passed as parameter if it varies
    capped_errors = [];
    if CAP_ERROR
        max_possible_error = sqrt(2) * area_size;
        capped_errors = metric.capErrorValues(all_errors, max_possible_error, INCLUDE_CAPPED);
        all_errors = capped_errors.values;
    end

    %% Select which metric to calculate and plot
    percentiles = struct('lower', [], 'upper', []); % Initialize percentiles for error band
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
        case 'band'
            percentiles = metric.cal_Percentiles(all_errors, BAND_PERCENTILES);
            plot_data = percentiles.val;
            metric_label = 'Median with Error Band';
        otherwise
            METRIC_TO_PLOT = 'p50';
            plot_data = metric.cal_Percentiles(all_errors).val;
            metric_label = 'Median';
    end
    
    %% Create display names and legends for all methods and plot the results
    legend_name = cell(1, num_legend);
    for i = 1:nvar_doa
        switch DOA_MODE
            case 'sweep'
                modeString = [DOA_MODE, ' ', num2str(DOA_RESOLUTION), '\circ res)'];
            case 'opt'
                modeString = [DOA_MODE, ' ', num2str(OPT_GRID_DENSITY), ' grid)'];
            otherwise
                modeString = DOA_MODE;
        end
        legend_name{i} = [strrep(doa_est_methods(i).name, '_', ' '), ' DoA triage (2 RXs ', modeString];
    end
    
    for ml_idx = 1:nvar_mlpos
        legend_name{nvar_doa+ml_idx} = ['MLpos ' num2str(RX_NUM(ml_idx)) ' RXs (triage initial)'];
        legend_name{nvar_doa+nvar_mlpos+ml_idx} = ['MLpos ' num2str(RX_NUM(ml_idx)) ' RXs (centroid initial)'];
    end
    
    rx_type = {'fixed', 'randomised'};
    cap_error = {'full', 'capped'};
    excluded = {' excluded', ''};
    annotStrings = {
        ['RX Type: ', rx_type{1 + RX_RANDOMISED}], ...
        ['ULA elements: ', num2str(ELEMENT_NUM)], ...
        ['Time instances: ', num2str(TIME_INST_NUM)], ...
        ['Error Metric: ', metric_label, ' (', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')']
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
    
    %% Plot time statistics for each algorithm 
    % Initialize data structures for each algorithm type
    ml_triage = struct('rx', [], 'time', [], 'std', []);
    ml_centroid = struct('rx', [], 'time', [], 'std', []);

    % Extract ML with triage initial data
    if isfield(timer_stats, 'triage_initial')
        rx_fields = fieldnames(timer_stats.triage_initial);
        for i = 1:length(rx_fields)
            rx_field = rx_fields{i};
            % Extract the RX count from "ML%dRXs" format
            rx_match = regexp(rx_field, 'ML(\d+)RXs', 'tokens');
            if ~isempty(rx_match)
                rx_count = str2double(rx_match{1}{1});
                ml_triage.rx(end+1) = rx_count;
                ml_triage.time(end+1) = timer_stats.triage_initial.(rx_field).mean;
                ml_triage.std(end+1) = timer_stats.triage_initial.(rx_field).std;
            end
        end
    end

    % Extract ML with centroid initial data
    if isfield(timer_stats, 'centroid_initial')
        rx_fields = fieldnames(timer_stats.centroid_initial);
        for i = 1:length(rx_fields)
            rx_field = rx_fields{i};
            % Extract the RX count from "ML%dRXs" format
            rx_match = regexp(rx_field, 'ML(\d+)RXs', 'tokens');
            if ~isempty(rx_match)
                rx_count = str2double(rx_match{1}{1});
                ml_centroid.rx(end+1) = rx_count;
                ml_centroid.time(end+1) = timer_stats.centroid_initial.(rx_field).mean;
                ml_centroid.std(end+1) = timer_stats.centroid_initial.(rx_field).std;
            end
        end
    end

    % Sort data by RX count
    [ml_triage.rx, idx] = sort(ml_triage.rx);
    ml_triage.time = ml_triage.time(idx);
    ml_triage.std = ml_triage.std(idx);

    [ml_centroid.rx, idx] = sort(ml_centroid.rx);
    ml_centroid.time = ml_centroid.time(idx);
    ml_centroid.std = ml_centroid.std(idx);

    % Create x_data and y_data for plotting
    x_data = [ml_triage.rx', ml_centroid.rx', ml_triage.rx']; % Last column for reference line
    y_data = [ml_triage.time', ml_centroid.time', repmat(timer_stats.DoA_triage.mean, size(ml_triage.rx'))];

    % Create error bands
    band_lower = [ml_triage.time' - ml_triage.std', ml_centroid.time' - ml_centroid.std', repmat(timer_stats.DoA_triage.mean - timer_stats.DoA_triage.std, size(ml_triage.rx'))];
    band_upper = [ml_triage.time' + ml_triage.std', ml_centroid.time' + ml_centroid.std', repmat(timer_stats.DoA_triage.mean + timer_stats.DoA_triage.std, size(ml_triage.rx'))];

    % Ensure no negative values in bands
    band_lower(band_lower < 0) = 0;

    % Plot execution times using the plots method
    metric.plots(x_data, y_data, 'linear', ...
        'DisplayNames', {' Triage initial ML', 'Centroid initial ML', 'DoA triage (reference)'}, ...
        'ShowBands', [true, true, false], ...
        'BandLower', band_lower, ...
        'BandUpper', band_upper, ...
        'Title', ['Time complexity Vs RX Count (', num2str(timer_stats.DoA_triage.count), ' iterations)'], ...
        'XLabel', 'Number of Receivers (RX)', ...
        'YLabel', 'Execution Time (seconds)', ...
        'LegendLocation', 'northwest', ...
        'ShowAnnotation', true);
                            
    % Optional: Add x-axis grid lines at each RX count
    grid on;
    xticks(unique([ml_triage.rx, ml_centroid.rx]));
    
    %% Plot empirical PDF of errors for selected SNR values
    snr_indices = [1, ceil(nvar_snr/2), nvar_snr];
    snr_values = mean(SNR_dB, 2);
    pdf_title = ['Empirical PDF of Position Estimation Errors by Method (',...
                rx_type{1 + RX_RANDOMISED},' RXs, ', num2str(ITERATION), ...
                ' iterations, ', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')'];
    
    metric.plot_epdf( ...
        all_errors, legend_name, snr_indices, pdf_title, ...
        'SNRValues', snr_values, ...
        'CapInfo', capped_errors, ...
        'CompareInSubplot', COMPARE_EPDF_IN_SUBPLOT);
        
    if CAP_ERROR
        % Display capped error values if applicable
        for idx = 1:num_legend
            fprintf('Method %d (%s): %d/%d values capped (%.2f%%)\n', idx, legend_name{idx}, ...
                capped_errors.cnt_capped{idx}, capped_errors.cnt_total{idx}, capped_errors.percentage{idx});
        end
    end
end

