clear; clc; close all;

%% User Inputs and Configurations
ITERATION = 500;                    % Number of Monte Carlo iterations
RX_CASES = [2,3];               % Different numbers of receivers to test
OPT_GRID_DENSITY = 10;              % Define a coarse grid for initial guesses
SNR_dB_RANGE = -10:2:20;            % SNR range in dB
ABS_ANGLE_LIM = 60;                 % Absolute angle limit (degrees)
TIME_INST_NUM = 1;                  % Number of time instances
RESOLUTION = 0.1;                   % Angle resolution (degrees)
FIXED_TRANS_ENERGY = true;          % Use fixed transmission energy
ELEMENT_NUM = 4;                    % Number of ULA elements
TX_SAFETY_DISTANCE = 2;             % Minimum distance between TX and RX (meters)

% Physical constants and wavelength
c = 299792458;                      % Speed of light (m/s)
fc = 2.4e9;                         % Operating frequency (Hz)
lambda = c / fc;                    % Wavelength

% Simulation area
area_size = 100;
pos_tx = [50, 50];                  % Tx at center
n_param = length(SNR_dB_RANGE);     % Number of SNR points to test

%% Initialize classes and arrays
channel = ChannelModels();
map2d = Map2D();
l4c = Likelihood4Coordinates();
optimiser = gridOptimiser();
rmse_all_cases = zeros(n_param, length(RX_CASES)); % Storage for results from different RX cases
%% Loop through each RX case
for rx_case_idx = 1:length(RX_CASES)
    RX_NUM = RX_CASES(rx_case_idx);  % Current number of receivers
    
    % Set SNR for all receivers 
    SNR_dB = repmat(SNR_dB_RANGE', 1, RX_NUM);
    
    % Display progress
    fprintf('Running simulation for %d receivers...\n', RX_NUM);
    progressbar('reset', ITERATION*n_param); % Reset progress bar
    
    %% Signal and channel configurations
    avg_amp_gain = 1;                % Average gain of the channel
    P_t = ones(RX_NUM, 1);           % W - Transmit signal power
    sub_carrier = (1:RX_NUM)' * 1000;% subcarrier spacing by 1000Hz
    Fs = 2 * max(sub_carrier);       % sample frequency
    T = TIME_INST_NUM/Fs;            % period of transmission
    t = 0:1/Fs:(T-1/Fs);             % Time vector for the signal
    
    % Receive Antenna elements characteristics
    element_spacing = 0.5 * lambda;  % Element spacing (ULA)
    
    %% Pre-calculate required values outside loop
    % Generate original transmitted signal
    s_t = sqrt(P_t(1)) .* exp(1j * 2 * pi * sub_carrier(1) * t);
    % Calculate average energy of the signal
    avg_E = FIXED_TRANS_ENERGY * 1 + ~FIXED_TRANS_ENERGY * (avg_amp_gain^2 * P_t(1) * T * Fs);
    
    %% Input
    y_los = channel.LoS(s_t, avg_amp_gain);
    ula = ULA(lambda, ELEMENT_NUM, element_spacing);
    rmse_values_ml = zeros(n_param, 1);
    w = cell(RX_NUM, 1); % Received signal at each Rx
    %% === Monte Carlo iterations
    for itr = 1:ITERATION
        %% --- Generate random RX positions ensuring minimum distance from TX
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
        
        % Generate random true Angle of Arrival
        aoa_act = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1);
        angle_rx_tx_abs = zeros(RX_NUM, 1);
        for i = 1:RX_NUM
            % Calculate the absolute angle of the receiver to the transmitter
            angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
        end
        rot_abs = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
        
        %% === Loop through each SNR value
        for snr_idx = 1:n_param
            % Update progress bar
            progressbar('step');
            %% === Generate the received signal at each Rx
            for rx_idx = 1:RX_NUM
                % --- Generate signal received at Rx
                nPower = avg_E/db2pow(SNR_dB(snr_idx, rx_idx));
                y_ula = channel.applyULA(y_los, aoa_act(rx_idx), ELEMENT_NUM, element_spacing, lambda);
                y_awgn = channel.AWGN(y_ula, nPower);
                % --- append received signal to a centralized array for ML estimation
                w{rx_idx} = y_awgn;
            end
            
            %% === Direct ML estimation
            objective_to_maximize = @(coor) -l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, w, ELEMENT_NUM, nPower);
            [optCoord, ~] = optimiser.fmincon2D(objective_to_maximize, {}, [0, 0], [area_size, area_size], OPT_GRID_DENSITY);
            rmse_values_ml(snr_idx) = rmse_values_ml(snr_idx) + sqrt((pos_tx(1,1)-optCoord(1))^2 + (pos_tx(1,2)-optCoord(2))^2);
        end
    end
    
    % Average RMSE values
    rmse_values_ml = rmse_values_ml / ITERATION;
    
    % Store results for this RX case
    rmse_all_cases(:, rx_case_idx) = rmse_values_ml;
end

%% === Plotting the comparison between different numbers of receivers
figure('Name', 'RMSE Comparison by Number of Receivers');
line_styles = {'-', '--', ':', '-.'};
markers = {'o', 's', 'd', '^'};
% colors = {'b', 'r', 'g', 'm'};

for i = 1:length(RX_CASES)
    semilogy(SNR_dB_RANGE, rmse_all_cases(:, i), ...
        [line_styles{mod(i-1, length(line_styles))+1}, markers{mod(i-1, length(markers))+1}], ...
        'LineWidth', 2, ...
        'MarkerSize', 8, ...
        'DisplayName', ['ML with ', num2str(RX_CASES(i)), ' receivers']); 
    grid on; hold on;
end

title(['RMSE Comparison by Number of Receivers (', num2str(ITERATION), ' iterations)']);
legend();
xlabel('Signal to Noise Ratio (SNR) [dB]');
ylabel('Root Mean Square Error (RMSE) [m]');

% Add a textbox with simulation parameters
annotation('textbox', [0.15, 0.1, 0.3, 0.2], ...
    'String', {['ULA Elements: ', num2str(ELEMENT_NUM)], ...
               ['Grid Density: ', num2str(OPT_GRID_DENSITY), 'x', num2str(OPT_GRID_DENSITY)], ...
               ['Area: ', num2str(area_size), 'x', num2str(area_size), ' m']}, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');