% filepath: /d:/workspaces/matlab/localisation/performance/snr_mse_MLpos_1tx_2rx.m
clear; clc; close all;

%% User Inputs and Configurations
ITERATION = 1000; % Reduced from 5000 for faster execution
OPT_GRID_DENSITY = 10; % Define a coarse grid for initial guesses
ABS_ANGLE_LIM = 0:30:60; % Different angle limits to test (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 0.1;                   % Angle resolution (degrees)
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements
DOA_MODE = 'sweep';                % DoA estimation mode ('sweep' or 'opt')
TX_SAFETY_DISTANCE = 2;             % Minimum distance between TX and RX (meters)
SHOW_ERROR_BAND = false;  % Whether to show the 25-75 percentile band
% Physical constants and wavelength
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength
% Transmitter, receiver positions angles
area_size = 100;
pos_tx = [50, 50];                  % Tx at center

RX_NUM = 2;                         % Number of receivers
SNR_dB = repmat((-10:2:20)', 1, RX_NUM);       % SNR in dB
n_param = length(SNR_dB); % Number of SNR points to test
n_angle_cases = length(ABS_ANGLE_LIM); % Number of angle limit cases



progressbar('reset', ITERATION*n_param*n_angle_cases); % Reset progress bar

%% Signal and channel configurations
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(RX_NUM, 1);  % W - Transmit signal power
sub_carrier = (1:RX_NUM)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA

% Generate original transmitted signal
s_t = sqrt(P_t(1)) .* exp(1j * 2 * pi * sub_carrier(1) * t);
% Calculate average energy of the signal
avg_E = FIXED_TRANS_ENERGY * 1 + ~FIXED_TRANS_ENERGY * (avg_amp_gain^2 * P_t(1) * T * Fs);

% Define the DoA estimation method
doa_est_method = 'MUSIC';
extra_args = {1};

% Add parameter to select which percentile to plot
PERCENTILE_TO_PLOT = 50; % Options: 25, 50 (median), 75

%% Replace mse_coor_val with all_errors for storing individual values
% Initialize storage for all individual errors
all_errors = cell(n_param, n_angle_cases); % Use cell array for variable-sized collections
for i=1:n_param
    for j=1:n_angle_cases
        all_errors{i,j} = zeros(ITERATION, 1); % Pre-allocate for each SNR-angle combination
    end
end

%% Initialise classes and arrays
channel = ChannelModels();
map2d = Map2D();
l4c = Likelihood4Coordinates();
optimiser = gridOptimiser();
y_los = channel.LoS(s_t, avg_amp_gain);
ula = ULA(lambda, ELEMENT_NUM, element_spacing);
w = cell(RX_NUM, 1); % Received signal at each Rx vectorised to cell array
% mse_coor_val = zeros(n_param, n_angle_cases); % Initialize storage for all angle limit cases
all_errors = cell(n_param, n_angle_cases); % Use cell array for variable-sized collections
for i=1:n_param
    for j=1:n_angle_cases
        all_errors{i,j} = zeros(ITERATION, 1); % Pre-allocate for each SNR-angle combination
    end
end
%% Loop through each angle limit case
for angle_idx = 1:n_angle_cases
    current_angle_limit = ABS_ANGLE_LIM(angle_idx);
    fprintf('Running simulation for angle limit: %d degrees...\n', current_angle_limit);
    
    %% === Monte Carlo iterations
    for itr = 1:ITERATION
        %% --- Location and AoA Refresh for each iteration
        % Generate random RX positions ensuring minimum distance from TX
        pos_rx = zeros(RX_NUM, 2);
        for i = 1:RX_NUM
            valid_position = false;
            while ~valid_position
                % Generate random position
                pos_rx(i,:) = area_size * rand(1, 2);

                % Check if it's far enough from TX
                if sqrt(sum((pos_tx - pos_rx(i,:)).^2)) >= TX_SAFETY_DISTANCE
                    valid_position = true;
                end
            end
        end
        
        % Generate random true Angle of Arrival within the current angle limit
        aoa_act = -current_angle_limit + RESOLUTION * randi([0, 2*current_angle_limit/RESOLUTION], RX_NUM, 1);
        angle_rx_tx_abs = zeros(RX_NUM, 1);
        for i = 1:RX_NUM
            % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
            angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
        end
        rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
        
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
                w{rx_idx} = y_awgn;
            end
            
            % --- DoA Estimation Algorithm at each RX
            aoa_rel_est = zeros(RX_NUM, 1);
            rays_abs = cell(RX_NUM, 1);
            
            for rx_idx = 1:RX_NUM
                estimator = DoAEstimator(ula, sweeping_angle, aoa_act(rx_idx), DOA_MODE, OPT_GRID_DENSITY);
                aoa_rel_est(rx_idx) = estimator.(doa_est_method)(w{rx_idx}, extra_args{:}).aoa_est;
                rays_abs{rx_idx} = map2d.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx));
            end
            
            % --- Calculate the aoa intersection point and the RMSE
            aoa_intersect = map2d.calDoAIntersect(rays_abs{1}, rays_abs{2});
            % Calculate error distance
            error_distance = sqrt((pos_tx(1,1)-aoa_intersect.x)^2 + (pos_tx(1,2)-aoa_intersect.y)^2);

            % Cap the error at a reasonable maximum (diagonal of area)
            max_possible_error = sqrt(2) * area_size; 
            if error_distance > max_possible_error || isnan(error_distance) || isinf(error_distance)
                error_distance = max_possible_error;
            end

            % mse_coor_val(snr_idx, angle_idx) = mse_coor_val(snr_idx, angle_idx) + error_distance^2;
            all_errors{snr_idx, angle_idx}(itr) = error_distance;
        end
    end
end
% Average the RMSE values for this angle limit case
% mse_coor_val = mse_coor_val / ITERATION;
% rmse_coor_val = sqrt(mse_coor_val);
percentile25 = zeros(n_param, n_angle_cases);
percentile50 = zeros(n_param, n_angle_cases);
percentile75 = zeros(n_param, n_angle_cases);
rmse_values = zeros(n_param, n_angle_cases);  % Keep RMSE for comparison

for i=1:n_param
    for j=1:n_angle_cases
        percentile25(i,j) = prctile(all_errors{i,j}, 25);
        percentile50(i,j) = prctile(all_errors{i,j}, 50);  % This is the median
        percentile75(i,j) = prctile(all_errors{i,j}, 75);
        rmse_values(i,j) = sqrt(mean(all_errors{i,j}.^2));  % Calculate RMSE too
    end
end

%% === Plotting
% Select which percentile to plot based on parameter
switch PERCENTILE_TO_PLOT
    case 25
        plot_data = percentile25;
        percentile_label = '25th Percentile';
    case 75
        plot_data = percentile75;
        percentile_label = '75th Percentile';
    otherwise
        plot_data = percentile50;  % Default to median
        percentile_label = 'Median (50th Percentile)';
end

figure('Name', ' Error Comparison by Angle Limit');

% Define line styles, markers and colors for different angle limits
line_styles = {'-', '--', ':', '-.', '-', '--'};
markers = {'o', 's', 'd', '^', 'v', 'p'};

% Plot selected percentile for each angle limit
for i = 1:n_angle_cases
    h_line = semilogy(mean(SNR_dB, 2), plot_data(:, i), ...
        [line_styles{mod(i-1, length(line_styles))+1}, markers{mod(i-1, length(markers))+1}], ...
        'LineWidth', 2, ...
        'MarkerSize', 6, ...
        'DisplayName', ['AoA Limit: \pm', num2str(ABS_ANGLE_LIM(i)), '\circ']); 
    grid on; hold on;
    
    % Optionally add error band (25th-75th percentile)
    if SHOW_ERROR_BAND && PERCENTILE_TO_PLOT == 50
        x = mean(SNR_dB, 2);
        % Get the color of the line that was just plotted
        line_color = get(h_line, 'Color');
        % Create shaded area between 25th and 75th percentiles
        fill([x; flipud(x)], ...
             [percentile25(:, i); flipud(percentile75(:, i))], ...
             line_color, ... % Use the same color as the line
             'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'DisplayName', ['Confidence Band \pm', num2str(ABS_ANGLE_LIM(i)), '\circ']);
    end
end

title(['Error Comparison by AoA Limit (', num2str(ITERATION), ' iterations)']);
legend('Location', 'northeast');
xlabel('Signal to Noise Ratio (SNR) [dB]');
ylabel([percentile_label, ' Error [m]']);

% Add a textbox with simulation parameters
annotation('textbox', [0.15, 0.1, 0.3, 0.2], ...
    'String', {['Method: ', strrep(doa_est_method, '_', ' ')], ...
               ['ULA Elements: ', num2str(ELEMENT_NUM)], ...
               ['Resolution: ', num2str(RESOLUTION), '\circ'], ...
               ['Error Metric: ', percentile_label]}, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');
