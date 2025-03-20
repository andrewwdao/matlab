clear; clc; close all;
%% User Inputs and Configurations
load('data/20250320_171313_aoa_mlpos_snr.mat');

%% === Prepare data and plotting
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
