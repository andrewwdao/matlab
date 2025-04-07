clear; clc; close all;

%% Select which data file to load
[filename, pathname] = uigetfile('data/*.mat', 'Select a saved simulation result');
if isequal(filename, 0)
    error('No file selected');
end

% Load the data file
load(fullfile(pathname, filename));

% Call the shared plotting function
plot_aoa_mlpos(all_errors, SNR_dB, metric, doa_est_methods, ...
    NUM_RX_DOA, NUM_RX_ML, DOA_MODE, DOA_RESOLUTION, OPT_GRID_DENSITY, ...
    ELEMENT_NUM, TIME_INST_NUM, RANDOMISE_RX, CAP_ERROR, INCLUDE_CAPPED, ...
    COMPARE_EPDF_IN_SUBPLOT, ITERATION, METRIC_TO_PLOT, BAND_PERCENTILES);
