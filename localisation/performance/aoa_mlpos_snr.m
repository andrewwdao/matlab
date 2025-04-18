%#ok<*UNRCH,*NASGU> % Suppress warnings for unreachable code and unused variables
clear; clc; close all;

%% User Inputs and Configurations
RUN_MODE = 'test';                  % Options: 'test', 'plot' or 'save'
METRIC_TO_PLOT = 'rmse';            % Options: 'rmse', 'p25', 'p50' (median), 'p75', 'band'
BAND_PERCENTILES = [25, 50, 75];    % Percentiles for error band if METRIC_TO_PLOT is 'band'
SHOW_ERROR_BAND = false;            % Whether to show the 25-75 percentile band
CAP_ERROR = false;                  % Cap error values at the maximum theoretical value
INCLUDE_CAPPED = true;              % Include capped values in the output errors, only valid if CAP_ERROR is true
COMPARE_EPDF_IN_SUBPLOT = false;    % Compare empirical PDFs in subplots
PLOT_METRICS = true;                % Plot metrics after simulation
PLOT_TIME_STATS = true;            % Plot time statistics after simulation
PLOT_EPDF = true;                   % Plot empirical PDF of errors after simulation

FLAG_PLOT = strcmp(RUN_MODE, 'test') || strcmp(RUN_MODE, 'plot'); % Flag for plotting
if strcmp(RUN_MODE, 'plot')
    % Get all variables that exist after line 11
    original_vars = who();
    
    % Store their values in a structure
    original_values = struct();
    for i = 1:length(original_vars)
        var_name = original_vars{i};
        original_values.(var_name) = evalin('base', var_name);
    end
    
    % Select which data file to load
    [filename, pathname] = uigetfile('data/*.mat', 'Select a saved simulation result');
    if isequal(filename, 0)
        error('No file selected');
    end
    
    % DON'T clear variables - just load the file directly
    % This automatically brings all file variables into the workspace
    load(fullfile(pathname, filename));
    
    % Now restore original configuration variables from original_values
    % This will override any variables from the file with the same names
    original_config_vars = fieldnames(original_values);
    for i = 1:length(original_config_vars)
        var_name = original_config_vars{i};
        assignin('base', var_name, original_values.(var_name));
    end
else
    % Set ITERATION based on run mode
    if strcmp(RUN_MODE, 'test')
        ITERATION = 1;  % For quick testing
        SAVE_METRICS = false;
    else  % 'full' mode
        ITERATION = 10000;  % For full simulation
        SAVE_METRICS = true;
    end

    TX_RANDOMISED = true;              % Randomise TX positions
    RX_RANDOMISED = true;              % Randomise RX positions and AoA
    TX_NUM = 1;                         % Number of transmitters
    RX_NUM = 3:7:24;                 % Additional receiver counts for ML optimization
    nvar_rx = length(RX_NUM);     % Number of variants for ML optimization
    ELEMENT_NUM = 4;                    % Number of ULA elements
    DOA_MODE = 'sweep';                 % DoA estimation mode ('sweep' or 'opt')
    DOA_RESOLUTION = 1;                 % Angle resolution (degrees)
    ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
    TIME_INST_NUM = 1;                  % Number of time instances
    FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
    OPT_GRID_DENSITY = 5;              % Define a coarse grid for initial guesses
    SAFETY_DISTANCE = 2;                % Minimum distance between TX and RX (meters)
    area_size = 100;                    % Area size for RX positions
    

    %% Initialise signal and channel configurations
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

    %% Initialise classes
    channel = ChannelModels();
    l4c = Likelihood4Coordinates();
    optimiser = Optimisers();
    algo = Algorithms(l4c, optimiser);
    map2d = Map2D([10,10], [90, 90], max(RX_NUM));
    metric = Metric();

    %% Initialise methods to test for performance
    methods_doa_est = struct(...
        'name', {'BF'}, ...            % estimator methods
        'extra_args', {{}} ...         % extra args required for specific type of estimator
    );
    nvar_doa_est = numel(methods_doa_est); % Get number of methods from the struct array
    % Define methods with merge configuration
    method_pos_est = struct(...
        'name', {'Centroid Direct ', 'Centroid MLpos ', 'Triage Direct ', 'Triage MLpos '}, ...
        'merge_rx', {false, false, true, false} ...  % Only merge Triage Direct across RX counts
    );
    nvar_pos_est = length(method_pos_est);
    legend4metric_num = nvar_rx*nvar_pos_est - sum([method_pos_est.merge_rx])*(nvar_rx-1); % Adjust legend count
    
    %% Initialise Transmitter, receiver positions and angles
    % Transmitters
    pos_tx = map2d.genTXPos(area_size, TX_NUM, TX_RANDOMISED);
    [s_t, e_avg] = channel.generateNuisanceSignal(fc, P_t, T, t, TIME_INST_NUM, FIXED_TRANS_ENERGY); % Generate nuisance transmitted signal with random phase
    % Initialise DoA estimator
    ula = ULA(lambda, ELEMENT_NUM, element_spacing);    % Create Uniform Linear Array object
    estimator = DoAEstimator(ula, sweeping_angle, 0, DOA_MODE, OPT_GRID_DENSITY);
    doa_estimator = @(sig) estimator.(methods_doa_est(1).name)(sig, methods_doa_est(1).extra_args{:});

    fprintf('Running Monte Carlo simulation with %d iterations...\n', ITERATION);
    progressbar('reset', ITERATION*nvar_snr+ ITERATION*nvar_rx*2*nvar_snr); % Reset progress bar
    %% Monte Carlo iterations
    % Preallocate error arrays for each method and SNR value
    errors_all = arrayfun(@(~) zeros(ITERATION, 1), ones(nvar_snr, legend4metric_num), 'UniformOutput', false);
    for itr = 1:ITERATION
        % --- Generate receivers and the received signal for the maximum number of receivers
        [pos_rx, aoa_act, rot_abs] = map2d.genRXPos(area_size, pos_tx, max(RX_NUM), RX_RANDOMISED, SAFETY_DISTANCE, ABS_ANGLE_LIM, DOA_RESOLUTION);
        [nPower, y_centralised] = channel.generateReceivedSignal(s_t, pos_tx, pos_rx, aoa_act, e_avg, SNR_dB, L_d0, d0, alpha, ELEMENT_NUM, element_spacing, lambda);
        
        %% --- Algorithms for position estimation
        for idx_rx = 1:nvar_rx
            % Create algorithm names with RX count
            timer_centroid = sprintf('centroid/MLpos%dRXs', RX_NUM(idx_rx));
            timer_triage = sprintf('triage/MLpos%dRXs', RX_NUM(idx_rx));
            % Loop through each SNR value
            for idx_snr=1:nvar_snr
                pos_rx_active = pos_rx(1:RX_NUM(idx_rx),:);
                y_received_active = y_centralised(idx_snr, 1:RX_NUM(idx_rx),:);
                progressbar('step'); % Update progress bar
                progressbar('starttimer', timer_centroid); % Time the algorithm
                [~, errors_all{idx_snr, 0*nvar_rx+idx_rx}(itr), ...
                 ~, errors_all{idx_snr, 1*nvar_rx+idx_rx}(itr), ~] = algo.MLOpt4mCentroid(...
                    pos_rx_active, rot_abs, y_received_active, ...
                    ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                    doa_estimator, pos_tx...
                );
                progressbar('stoptimer', timer_centroid);
                progressbar('step'); % Update progress bar
                progressbar('starttimer', timer_triage); % Time the algorithm
                [~, errors_all{idx_snr, 2*nvar_rx+idx_rx}(itr), ...
                 ~, errors_all{idx_snr, 3*nvar_rx+idx_rx}(itr), ~] = algo.MLOpt4mDoAtriage(...
                    pos_rx_active, rot_abs, y_received_active, ...
                    ELEMENT_NUM, nPower, [0, 0], [area_size, area_size],...
                    doa_estimator, pos_tx...
                );
                progressbar('stoptimer', timer_triage);
                % progressbar('step'); % Update progress bar
                % [~, errors_all{idx_snr, nvar_doa_est+2*nvar_rx+idx_rx}(itr), ~] = algo.MLOptwGrid(...
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
            fprintf('Creating directory: %s\n', OUTPUT_PATH);
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
    capped_errors = [];
    if CAP_ERROR
        max_possible_error = sqrt(2) * area_size;
        capped_errors = metric.capErrorValues(errors_all, max_possible_error, INCLUDE_CAPPED);
        errors_all = capped_errors.values;
        % Display capped error values
        for idx = 1:legend4metric_num
            fprintf('Method %d (%s): %d/%d values capped (%.2f%%)\n', idx, legend4metric_name{idx}, ...
                capped_errors.cnt_capped{idx}, capped_errors.cnt_total{idx}, capped_errors.percentage{idx});
        end
    end

    %% --- Merge errors for identical methods
    % Initialize a counter to track the total number of columns needed after merging
    total_cols = 0;
    for idx_method = 1:nvar_pos_est
        if method_pos_est(idx_method).merge_rx
            total_cols = total_cols + 1;  % One column for merged method
        else
            total_cols = total_cols + nvar_rx;  % One column per RX count
        end
    end

    % Initialize properly sized merged errors array
    errors_merged = cell(size(errors_all, 1), total_cols);

    for idx_snr = 1:size(errors_all, 1)
        merged_col = 1; % Track position in the new merged array
        
        for idx_method = 1:nvar_pos_est
            if method_pos_est(idx_method).merge_rx
                % For methods to merge, average across all RX counts
                merged_data = zeros(size(errors_all{1}));
                
                % Calculate average across all RX counts for this method
                for idx_rx = 1:nvar_rx
                    orig_idx = (idx_method-1)*nvar_rx + idx_rx;
                    merged_data = merged_data + errors_all{idx_snr, orig_idx};
                end
                
                % Store the merged data in the proper position
                errors_merged{idx_snr, merged_col} = merged_data / nvar_rx;
                merged_col = merged_col + 1;
            else
                % For other methods, keep separate entries for each RX count
                for idx_rx = 1:nvar_rx
                    orig_idx = (idx_method-1)*nvar_rx + idx_rx;
                    errors_merged{idx_snr, merged_col} = errors_all{idx_snr, orig_idx};
                    merged_col = merged_col + 1;
                end
            end
        end
    end

    %% --- Calculate selected metric
    percentiles = struct('lower', [], 'upper', []); % Initialize percentiles for error band
    switch METRIC_TO_PLOT
        case 'rmse'
            metric_plot_data = metric.cal_RMSE(errors_merged);
            metric_label = 'RMSE';
        case 'p25'
            metric_plot_data = metric.cal_Percentiles(errors_merged, 25).val;
            metric_label = '25th Percentile';
        case 'p75'
            metric_plot_data = metric.cal_Percentiles(errors_merged, 75).val;
            metric_label = '75th Percentile';
        case 'band'
            percentiles = metric.cal_Percentiles(errors_merged, BAND_PERCENTILES);
            metric_plot_data = percentiles.val;
            metric_label = 'Median with Error Band';
        otherwise
            METRIC_TO_PLOT = 'p50';
            metric_plot_data = metric.cal_Percentiles(errors_merged).val;
            metric_label = 'Median';
    end
    
    %% --- Create annotation strings and legend for the plots
    % Annotation strings
    gen_type = {'fixed', 'randomised'};
    cap_error = {'full', 'capped'};
    excluded = {' excluded', ''};
    switch DOA_MODE
        case 'sweep'
            modeString = [DOA_MODE ' ' num2str(DOA_RESOLUTION) '\circ res'];
        case 'opt'
            modeString = [DOA_MODE ' ' num2str(OPT_GRID_DENSITY) ' grid'];
        otherwise
            modeString = DOA_MODE;
    end
    annotStrings = {
        [gen_type{1 + TX_RANDOMISED}, ' TX, ', gen_type{1 + RX_RANDOMISED}, ' RX'], ...
        ['ULA elements: ', num2str(ELEMENT_NUM)], ...
        ['Time instances: ', num2str(TIME_INST_NUM)], ...
        ['DoA method: ', modeString], ...
        ['Error Metric: ', metric_label, ' (', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')']
    };

    %% --- Plot the error metric for each algorithm
    if PLOT_METRICS
        % Create display names and legends for all methods and plot the results
        legend4metric_name = cell(1, legend4metric_num);
        idx_legend = 1; % Track the current position in the legend array

        % Initialize the plot data for the merged metrics
        for idx_method = 1:nvar_pos_est
            if method_pos_est(idx_method).merge_rx
                % For methods to merge, add a single legend entry
                legend4metric_name{idx_legend} = [method_pos_est(idx_method).name, '(x', num2str(nvar_rx), ' itr)'];
                idx_legend = idx_legend + 1;
            else
                % For other methods, keep separate entries for each RX count
                for idx_rx = 1:nvar_rx
                    orig_idx = (idx_method-1)*nvar_rx + idx_rx;
                    legend4metric_name{idx_legend} = [method_pos_est(idx_method).name, num2str(RX_NUM(idx_rx)), ' RXs'];
                    idx_legend = idx_legend + 1;
                end
            end
        end

        % Plot with the merged data
        metric.plots(mean(SNR_dB, 2), metric_plot_data, 'semilogy', ...
            'DisplayNames', legend4metric_name, ...
            'ShowBands', strcmp(METRIC_TO_PLOT, 'band') * ones(1, legend4metric_num), ...
            'BandLower', percentiles.lower, ...
            'BandUpper', percentiles.upper, ...
            'Title', ['Error Metric by estimation method (', num2str(ITERATION), ' iterations)'], ...
            'YLabel', [metric_label, ' Error [m]'], ...
            'ShowAnnotation', true, ...
            'AnnotationStrings', annotStrings);
    end
    %% --- Plot time statistics for each algorithm
    if PLOT_TIME_STATS
        % Extract unique parent algorithm names from the timer data
        unique_parents = unique({timer_stats.Parent});
        % Remove empty strings if any top-level timers exist without children defined this way
        unique_parents = unique_parents(~cellfun('isempty', unique_parents)); 

        % Prepare data structures for plotting
        max_entries = length(timer_stats); % Maximum possible unique combinations
        timer_stats_plot_data = repmat(struct('Parent', '', 'Type', '', 'RX', [], 'Mean', [], 'Std', []), 1, max_entries);
        timer_stats_cnt = 0; % Track actual used entries
        itr_cnt = 0; % To store the number of iterations

        % --- Aggregate data from timer_stats ---
        for i = 1:length(timer_stats)
            entry = timer_stats(i);
            
            % Try to parse Type (DoA/ML) and RX count from the Name field
            rx_match = regexp(entry.Name, '^(Direct|MLpos)(\d+)RXs$', 'tokens', 'once'); % Match DoA or ML
            
            if ~isempty(rx_match) && isfield(entry, 'Parent') && ~isempty(entry.Parent)
                timer_type = rx_match{1};
                rx_count = str2double(rx_match{2});
                parent_name = entry.Parent;

                % Find if this parent/type combo already exists in timer_stats_plot_data
                found = false;
                for k = 1:timer_stats_cnt
                    if strcmp(timer_stats_plot_data(k).Parent, parent_name) && strcmp(timer_stats_plot_data(k).Type, timer_type)
                        timer_stats_plot_data(k).RX(end+1) = rx_count;
                        timer_stats_plot_data(k).Mean(end+1) = entry.Mean;
                        timer_stats_plot_data(k).Std(end+1) = entry.Std;
                        found = true;
                        break;
                    end
                end
                
                % If not found, add a new entry to timer_stats_plot_data
                if ~found
                    timer_stats_cnt = timer_stats_cnt + 1;
                    timer_stats_plot_data(timer_stats_cnt) = struct(...
                        'Parent', parent_name, ...
                        'Type', timer_type, ...
                        'RX', rx_count, ...
                        'Mean', entry.Mean, ...
                        'Std', entry.Std);
                end

                if itr_cnt == 0 && isfield(entry, 'Count')
                    itr_cnt = entry.Count;
                end
            end
        end

        % Trim unused entries
        timer_stats_plot_data = timer_stats_plot_data(1:timer_stats_cnt);

        % --- Prepare final arrays for plotting from aggregated timer_stats_plot_data ---
        % First calculate the sizes needed for each array
        x_data = cell(1, timer_stats_cnt);
        y_data = cell(1, timer_stats_cnt);
        lower_data = cell(1, timer_stats_cnt);
        upper_data = cell(1, timer_stats_cnt);
        final_legend_names = cell(1, timer_stats_cnt);

        for k = 1:length(timer_stats_plot_data)
            % Sort data by RX count for plotting
            [sorted_rx, sort_idx] = sort(timer_stats_plot_data(k).RX);
            sorted_mean = timer_stats_plot_data(k).Mean(sort_idx);
            sorted_std = timer_stats_plot_data(k).Std(sort_idx);
            
            % Store in cell arrays
            x_data{k} = sorted_rx';
            y_data{k} = sorted_mean';
            lower_data{k} = max(sorted_mean' - sorted_std', 0); % Ensure non-negative
            upper_data{k} = sorted_mean' + sorted_std';
            
            % Generate legend name
            parent_display = strrep(timer_stats_plot_data(k).Parent, '_', ' ');
            parent_display = [upper(parent_display(1)), parent_display(2:end)]; % Capitalize
            final_legend_names{k} = sprintf('%s %s', parent_display, timer_stats_plot_data(k).Type);
        end

        % Convert cell arrays to matrices for plotting
        x_plot = cell2mat(x_data);
        y_plot = cell2mat(y_data);
        band_lower_plot = cell2mat(lower_data);
        band_upper_plot = cell2mat(upper_data);

        % Check if we have data to plot
        if ~isempty(x_plot)
            % Plot time statistics
            metric.plots(x_plot, y_plot, 'linear', ...
                'DisplayNames', final_legend_names, ...
                'ShowBands', true(1, size(x_plot, 2)), ... % Show bands for all lines
                'BandLower', band_lower_plot, ...
                'BandUpper', band_upper_plot, ...
                'Title', ['Time complexity Vs RX Count (', num2str(itr_cnt), ' iterations)'], ...
                'XLabel', 'Number of Receivers (RX)', ...
                'YLabel', 'Execution Time (seconds)', ...
                'LegendLocation', 'north', ... 
                'ShowAnnotation', true, ...
                'AnnotationPosition', [0.14, 0.7, 0.3, 0.2], ... 
                'AnnotationStrings', annotStrings); 
        else
            % Updated error message for clarity
            fprintf('No timer data found matching the expected format (e.g., Parent/DoA3RXs, Parent/ML10RXs) for plotting time complexity vs RX count.\n');
        end
    end
    %% --- Plot empirical PDF of errors for selected SNR values
    if PLOT_EPDF
        snr_indices = [1, ceil(nvar_snr/2), nvar_snr];
        snr_values = mean(SNR_dB, 2);
        pdf_title = ['Empirical PDF of Position Estimation Errors by Method (',...
                    gen_type{1 + RX_RANDOMISED},' RXs, ', num2str(ITERATION), ...
                    ' iterations, ', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')'];
        
        metric.plot_epdf( ...
            errors_merged, legend4metric_name, snr_indices, pdf_title, ...
            'SNRValues', snr_values, ...
            'CapInfo', capped_errors, ...
            'CompareInSubplot', COMPARE_EPDF_IN_SUBPLOT);
    end
end
