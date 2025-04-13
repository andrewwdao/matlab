%#ok<*UNRCH,*NASGU> % Suppress warnings for unreachable code and unused variables
clear; clc; close all;

%% User Inputs and Configurations
RUN_MODE = 'plot';                  % Options: 'test', 'plot' or 'save'
METRIC_TO_PLOT = 'rmse';            % Options: 'rmse', 'p25', 'p50' (median), 'p75', 'band'
BAND_PERCENTILES = [25, 50, 75];    % Percentiles for error band if METRIC_TO_PLOT is 'band'
SHOW_ERROR_BAND = false;            % Whether to show the 25-75 percentile band
CAP_ERROR = false;                  % Cap error values at the maximum theoretical value
INCLUDE_CAPPED = true;              % Include capped values in the output errors, only valid if CAP_ERROR is true
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
        ITERATION = 10000;  % For full simulation
        SAVE_METRICS = true;
    end

    TX_RANDOMISED = false;              % Randomise TX positions
    RX_RANDOMISED = false;              % Randomise RX positions and AoA
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
    methods_doa_est = struct(...
        'name', {'BF'}, ...            % estimator methods
        'extra_args', {{}} ...         % extra args required for specific type of estimator
    );
    nvar_doa_est = numel(methods_doa_est); % Get number of methods from struct array
    method_pos_est = {
        'Triage ', 'MLpos Triage ', 'Centroid ', 'MLpos Centroid '
    };
    nvar_pos_est = length(method_pos_est); % Get number of methods from cell array
    legend4metric_num = nvar_mlpos*nvar_pos_est;         % Number of methods plus ML method

    % Initialize storage for all individual errors
    all_errors = cell(nvar_snr, legend4metric_num);             % Pre-allocate for variable-sized collections
    for i=1:nvar_snr
        for j=1:legend4metric_num
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
    doa_estimator = @(sig) estimator.(methods_doa_est(1).name)(sig, methods_doa_est(1).extra_args{:});

    fprintf('Running Monte Carlo simulation with %d iterations...\n', ITERATION);
    progressbar('reset', ITERATION*nvar_snr+ ITERATION*nvar_mlpos*2*nvar_snr); % Reset progress bar
    %% === Monte Carlo iterations
    for itr = 1:ITERATION
        % --- Generate receivers and the received signal for the maximum number of receivers
        [pos_rx, aoa_act, rot_abs] = map2d.genRXPos(area_size, pos_tx, max(RX_NUM), RX_RANDOMISED, SAFETY_DISTANCE, ABS_ANGLE_LIM, DOA_RESOLUTION);
        [nPower, y_centralised] = channel.generateReceivedSignal(s_t, pos_tx, pos_rx, aoa_act, e_avg, SNR_dB, L_d0, d0, alpha, ELEMENT_NUM, element_spacing, lambda);

        % %% --- DoA estimation only
        % for idx_snr=1:nvar_snr
        %     y_received = y_centralised(idx_snr, :);
        %     progressbar('step'); % Update progress bar
        %     % Time the DoAtriage algorithm
        %     progressbar('starttimer', 'DoA_triage');
        %     [~, all_errors{idx_snr, 1}(itr)] = algo.DoAtriage(...
        %         pos_rx, rot_abs, y_received, ...
        %         doa_estimator, pos_tx ...
        %     );
        %     progressbar('stoptimer', 'DoA_triage');
        % end
        
        %% --- ML optimization with additional receivers
        for ml_idx = 1:nvar_mlpos
            % Create algorithm names with RX count
            triage_timer_name = sprintf('triage/ML%dRXs', RX_NUM(ml_idx));
            centroid_timer_name = sprintf('centroid/ML%dRXs', RX_NUM(ml_idx));
            % Loop through each SNR value
            for idx_snr=1:nvar_snr
                pos_rx_active = pos_rx(1:RX_NUM(ml_idx),:);
                y_received_active = y_centralised(idx_snr, 1:RX_NUM(ml_idx),:);
                progressbar('step'); % Update progress bar
                progressbar('starttimer', triage_timer_name); % Time the algorithm
                [~, all_errors{idx_snr, ml_idx}(itr), ~, all_errors{idx_snr, nvar_mlpos+ ml_idx}(itr), ~] = algo.MLOpt4mDoAtriage(...
                    pos_rx_active, rot_abs, y_received_active, ...
                    ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                    doa_estimator, pos_tx...
                );
                progressbar('stoptimer', triage_timer_name);
                progressbar('step'); % Update progress bar
                progressbar('starttimer', centroid_timer_name); % Time the algorithm
                [~, all_errors{idx_snr, 2*nvar_mlpos+ml_idx}(itr), ~, all_errors{idx_snr, 3*nvar_mlpos+ml_idx}(itr), ~] = algo.MLOpt4mCentroid(...
                    pos_rx_active, rot_abs, y_received_active, ...
                    ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                    doa_estimator, pos_tx...
                );
                progressbar('stoptimer', centroid_timer_name);
                % progressbar('step'); % Update progress bar
                % [~, all_errors{idx_snr, nvar_doa_est+2*nvar_mlpos+ml_idx}(itr), ~] = algo.MLOptwGrid(...
                %     pos_rx_active, rot_abs, y_received_active, ...
                %     ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                %     OPT_GRID_DENSITY, pos_tx...
                % );
            end
        end
    end
% After all iterations, report the execution times
timer_stats = progressbar('reporttimers');
    %% --- Save metrics to file
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
    %% --- Cap errors at the maximum theoretical value
    area_size = 100; % This should be passed as parameter if it varies
    capped_errors = [];
    if CAP_ERROR
        max_possible_error = sqrt(2) * area_size;
        capped_errors = metric.capErrorValues(all_errors, max_possible_error, INCLUDE_CAPPED);
        all_errors = capped_errors.values;
    end

    %% --- Calculate selected metric
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
    
    %% --- Create annotation strings and legend for the plots
    % Annotation strings
    rx_type = {'fixed', 'randomised'};
    cap_error = {'full', 'capped'};
    excluded = {' excluded', ''};
    annotStrings = {
        ['RX Type: ', rx_type{1 + RX_RANDOMISED}], ...
        ['ULA elements: ', num2str(ELEMENT_NUM)], ...
        ['Time instances: ', num2str(TIME_INST_NUM)], ...
        ['Error Metric: ', metric_label, ' (', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')']
    };

    %% --- Plot the error metric for each algorithm 
    % Create display names and legends for all methods and plot the results
    legend4metric_name = cell(1, legend4metric_num);
    modeString = DOA_MODE;
    if strcmp(DOA_MODE, 'sweep')
        modeString = [modeString, ' ', num2str(DOA_RESOLUTION), '\circ res'];
    elseif strcmp(DOA_MODE, 'opt')
        modeString = [modeString, ' ', num2str(OPT_GRID_DENSITY), ' grid'];
    end
    
    % Create legend names for DoA methods
    % for i = 1:nvar_doa_est
    %     legend4metric_name{i} = [method_pos_est{1}, ' 2 RXs, ', modeString, ', ', strrep(methods_doa_est(i).name, '_', ' ')];
    % end
    
    % Create legend names for ML methods dynamically
    num_ml_methods = length(method_pos_est);
    for method_idx = 1:num_ml_methods
        offset = (method_idx-1)*nvar_mlpos;
        for ml_idx = 1:nvar_mlpos
            legend4metric_name{offset + ml_idx} = [method_pos_est{method_idx}, num2str(RX_NUM(ml_idx)) ' RXs'];
        end
    end

    % Plot the error metric
    metric.plots(mean(SNR_dB, 2), plot_data, 'semilogy', ...
        'DisplayNames', legend4metric_name, ...
        'ShowBands', strcmp(METRIC_TO_PLOT, 'band') * ones(1, legend4metric_num), ...
        'BandLower', percentiles.lower, ...
        'BandUpper', percentiles.upper, ...
        'Title', ['Error Metric by estimation method (', num2str(ITERATION), ' iterations)'], ...
        'YLabel', [metric_label, ' Error [m]'], ...
        'ShowAnnotation', true, ...
        'AnnotationStrings', annotStrings);
    
    %% --- Plot time statistics for each algorithm 
    % Extract data for both algorithm types
    ml_triage = extractAlgorithmData(timer_stats, 'triage');
    ml_centroid = extractAlgorithmData(timer_stats, 'centroid');

    % % Extract reference data once
    % ref_size = size(ml_triage.rx');
    % ref_mean = repmat(timer_stats.triage.DoA2RXs.mean, ref_size);
    % ref_std = repmat(timer_stats.triage.DoA2RXs.std, ref_size);

    % Prepare all data arrays
    x_data = [ml_triage.rx', ml_centroid.rx'];
    means = [ml_triage.mean', ml_centroid.mean'];
    stds = [ml_triage.std', ml_centroid.std'];

    % Create plotting data
    y_data = means;
    band_lower = max(means - stds, 0); % Prevents negative values
    band_upper = means + stds;

    % Plot time statistics
    metric.plots(x_data, y_data, 'linear', ...
        'DisplayNames', method_pos_est, ...
        'ShowBands', [true, true, true, true], ...
        'BandLower', band_lower, ...
        'BandUpper', band_upper, ...
        'Title', ['Time complexity Vs RX Count (', num2str(timer_stats.triage.ML3RXs.count), ' iterations)'], ...
        'XLabel', 'Number of Receivers (RX)', ...
        'YLabel', 'Execution Time (seconds)', ...
        'LegendLocation', 'north', ...
        'ShowAnnotation', true, ...
        'AnnotationPosition', [0.14, 0.7, 0.3, 0.2], ...
        'AnnotationStrings', annotStrings);
    
    %% --- Plot empirical PDF of errors for selected SNR values
    snr_indices = [1, ceil(nvar_snr/2), nvar_snr];
    snr_values = mean(SNR_dB, 2);
    pdf_title = ['Empirical PDF of Position Estimation Errors by Method (',...
                rx_type{1 + RX_RANDOMISED},' RXs, ', num2str(ITERATION), ...
                ' iterations, ', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')'];
    
    metric.plot_epdf( ...
        all_errors, legend4metric_name, snr_indices, pdf_title, ...
        'SNRValues', snr_values, ...
        'CapInfo', capped_errors, ...
        'CompareInSubplot', COMPARE_EPDF_IN_SUBPLOT);
        
    if CAP_ERROR
        % Display capped error values if applicable
        for idx = 1:legend4metric_num
            fprintf('Method %d (%s): %d/%d values capped (%.2f%%)\n', idx, legend4metric_name{idx}, ...
                capped_errors.cnt_capped{idx}, capped_errors.cnt_total{idx}, capped_errors.percentage{idx});
        end
    end
end

% Function to extract and sort algorithm data with DoA first then ML
function data = extractAlgorithmData(timer_stats, field_name)
    data = struct('rx', [], 'mean', [], 'std', [], 'method', [], 'field_name', []);
    
    if isfield(timer_stats, field_name)
        method_fields = fieldnames(timer_stats.(field_name));
        for i = 1:length(method_fields)
            current_field = method_fields{i};
            % Extract method type (DoA or ML) and RX count
            method_match = regexp(current_field, '^(DoA|ML)(\d+)RXs', 'tokens', 'once');
            if ~isempty(method_match)
                method_type = method_match{1};
                rx_count = str2double(method_match{2});
                data.rx(end+1) = rx_count;
                data.method{end+1} = method_type;
                data.field_name{end+1} = current_field;
                data.mean(end+1) = timer_stats.(field_name).(current_field).mean;
                data.std(end+1) = timer_stats.(field_name).(current_field).std;
            end
        end
        
        % Custom sorting: DoA first, then ML, then by RX count
        if ~isempty(data.rx)
            % Create a sorting value: DoA=0, ML=1 (for primary sort by method)
            method_values = zeros(size(data.method));
            for i = 1:length(data.method)
                if strcmp(data.method{i}, 'ML')
                    method_values(i) = 1;
                end
            end
            
            % Create a composite sorting key: method_value*1000 + rx_count
            % This ensures DoA comes before ML, and within each method type, sorted by RX count
            sort_keys = method_values*1000 + data.rx;
            [~, idx] = sort(sort_keys);
            
            % Reorder all fields based on the sorted indices
            data.rx = data.rx(idx);
            data.mean = data.mean(idx);
            data.std = data.std(idx);
            data.method = data.method(idx);
            data.field_name = data.field_name(idx);
        end
    end
end