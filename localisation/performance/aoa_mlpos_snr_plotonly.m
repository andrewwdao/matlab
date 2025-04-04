clear; clc; close all;
%% User Inputs and Configurations
load('data/20250329_150352_aoa_mlpos_snr.mat');

%% === Prepare data and plotting
% Cap errors at the maximum theoretical value
if CAP_ERROR
    max_possible_error = sqrt(2) * area_size;
    capped_errors = metric.capErrorValues(all_errors, max_possible_error, INCLUDE_CAPPED);
    all_errors = capped_errors.values;
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
            modeString = [DOA_MODE, ' ', num2str(DOA_RESOLUTION), '\circ res)'];
        case 'opt'
            modeString = [DOA_MODE, ' ', num2str(OPT_GRID_DENSITY), ' grid)'];
        otherwise
            modeString = DOA_MODE;
    end
    legend_name{i} = [strrep(doa_est_methods(i).name, '_', ' '), ' DoA triage (', num2str(NUM_RX_DOA), ' RXs ', modeString];
end
for ml_idx = 1:nvar_mlpos
    legend_name{nvar_doa+ml_idx} = ['MLpos ' num2str(NUM_RX_ML(ml_idx)) ' RXs (triage initial)'];
    legend_name{nvar_doa+nvar_mlpos+ml_idx} = ['MLpos ' num2str(NUM_RX_ML(ml_idx)) ' RXs (centroid initial)'];
end

rx_type = {'fixed', 'randomised'};
cap_error = {'full', 'capped'};
excluded = {' excluded', ''};
annotStrings = {
    ['RX Type: ', rx_type{1 + RANDOMISE_RX}], ...
    ['ULA elements: ', num2str(ELEMENT_NUM)], ...
    ['Time instances: ', num2str(TIME_INST_NUM)], ...
    ['Error Metric: ', metric_label, ' (', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')']
};

% Display capped error values
if CAP_ERROR
    for idx =1:num_legend
        fprintf('Method %d (%s): %d/%d values capped (%.2f%%)\n', idx, legend_name{idx}, capped_errors.cnt_capped{idx}, capped_errors.cnt_total{idx}, capped_errors.percentage{idx});
    end
end
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

%% Plot empirical PDF of errors for selected SNR values
% Select low, medium, and high SNR indices to plot
snr_indices = [1, ceil(nvar_snr/2), nvar_snr];

% Get SNR values for the selected indices
snr_values = mean(SNR_dB, 2);

% Add a figure title
pdf_title = ['Empirical PDF of Position Estimation Errors by Method (',rx_type{1 + RANDOMISE_RX},' RXs, ', num2str(ITERATION), ' iterations, ', cap_error{1+CAP_ERROR}, excluded{1+INCLUDE_CAPPED},')'];

% Call the enhanced plot_epdf method from the Metric class
if CAP_ERROR
    metric.plot_epdf(all_errors, legend_name, snr_indices, pdf_title, ...
                    'SNRValues', snr_values, ...
                    'CapInfo', capped_errors, ...
                    'ShowKDE', true, ...
                    'CompareMethod', 1, ... % Compare all methods to DoA intersection
                    'CompareInSubplot', COMPARE_EPDF_IN_SUBPLOT);
else
    metric.plot_epdf(all_errors, legend_name, snr_indices, pdf_title, ...
                    'SNRValues', snr_values, ...
                    'ShowKDE', true, ...
                    'CompareMethod', 1, ... % Compare all methods to DoA intersection
                    'CompareInSubplot', COMPARE_EPDF_IN_SUBPLOT);
end