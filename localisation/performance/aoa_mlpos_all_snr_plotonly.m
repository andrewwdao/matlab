clear; clc; close all;
%% User Inputs and Configurations
load('data/20250318_123605_snr_emetric_aoa_n_mlpos_plotonly.mat');
%% Initialise classes
channel = ChannelModels();
map2d = Map2D();
metric = Metric();
l4c = Likelihood4Coordinates();
optimiser = Optimisers();
%% SNR values to test
nvar_snr = length(SNR_dB);                   % Number of positions to test
nvar_mlpos = length(NUM_RX_ML);     % Number of variants for ML optimization

%% Signal and channel configurations
nvar_doa = numel(doa_est_methods);               % Automatically get number of methods from struct array
num_legend = nvar_doa + nvar_mlpos;         % Number of methods plus ML method

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
