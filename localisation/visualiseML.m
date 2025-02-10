clear; clc; close all;
%% visualize_L - Computes and plots L(x,y) = P(f1(x,y), f2(x,y)) in 3D.
%% === User inputs
RX_NUM = 2; % Number of receivers
SNR_dB = -10*ones(RX_NUM, 1); %dB
SHOW_LIMITS = true; % Show the detecting limits of the RXs (with known limitation)
ABS_ANGLE_LIM = 60; % Absolute angle limit in degree
TIME_INST_NUM = 1; % Number of time instances
RESOLUTION = 0.1; % Angle resolution in degreetranslate to -lim to lim (symmetric)
FIXED_TRANS_ENERGY = true; % Flag to use Average SNR over all time instances or SNR over ONE time instance
ELEMENT_NUM = 4; % Number of elements in the ULA
%% === Other configurations
% rs=rng(2007); % initialize the random number generator to a specific seed value
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
% --- Location Calculation
area_size = 100;   % 100x100 meter area
% aoa_true = -ABS_ANGLE_LIM + RESOLUTION * randi([0, 2*ABS_ANGLE_LIM/RESOLUTION], RX_NUM, 1); % true Angle of Arrival from RX to TX that will later be transformed to the absolute angle
% aoa_true = [-6.8; 45.6];
aoa_true = [-1.1; -5.4];
pos_tx = [50, 50]; % Transmitter position (x, y) in meters
% pos_rx = area_size*rand(RX_NUM, 2); % Random Receiver position (x, y) in meters
pos_rx = [15.98, 54.27; 67.64, 69.28];

angle_rx_tx_abs = zeros(RX_NUM, 1);
for i = 1:RX_NUM  % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
    angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
end
rot_abs = angle_rx_tx_abs - aoa_true; % Absolute rotation of the receiver in degrees

progressbar('reset', RX_NUM); % Reset progress bar
progressbar('displaymode', 'append'); % Reset progress bar
progressbar('minimalupdateinterval', 1); % Update progress bar every x seconds
avg_amp_gain = 1; % Average gain of the channel
P_t = ones(RX_NUM, 1);  % W - Transmit signal power
sub_carrier = (1:RX_NUM)' * 1000;  % subcarrier spacing by 1000Hz
Fs = 2 * max(sub_carrier);  % sample frequency
T = TIME_INST_NUM/Fs; % period of transmission
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
% --- Receive Antenna elements characteristics
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
sweeping_angle = -90:RESOLUTION:90; % Angle range for finding the AoA
%% === Loop through each RX
% aoa_rel_est = zeros(RX_NUM, 1);
% rays_abs = cell(RX_NUM, 1);
% estimator_coor = PosEstimator2D();
w = cell(RX_NUM, 1);
for rx_idx=1:RX_NUM
    progressbar('advance'); % Update progress bar
    % Generate base signal
    s_t = sqrt(P_t(RX_NUM)) .* exp(1j * 2 * pi * sub_carrier(RX_NUM) * t);
    % Calculate average energy of the signal
    if FIXED_TRANS_ENERGY == true
        % Average engery is fixed for whole transmission time
        % regardless of the number of time instances.
        avg_E = 1;
    else
        % Average engery is fixed for one time instance,
        % so the whole energy over all time instances need to be calculated.
        avg_E = avg_amp_gain^2 * P_t(RX_NUM) * T * Fs;
    end
    % Calculate noise parameters with the corresponding average energy and SNR
    nPower = avg_E/db2pow(SNR_dB(rx_idx));
    % Initialize channel model
    channel = ChannelModelAoA(aoa_true(rx_idx), lambda, ELEMENT_NUM, element_spacing);
    %% === Generate original signal received at Rx
    y_los = channel.LoS(s_t, avg_amp_gain);  % Received signal at the receiver
    y_ula = channel.applyULA(y_los);  % Apply ULA characteristics to the received signal
    y_awgn = channel.AWGN(y_ula, nPower);
    w{rx_idx} = y_awgn;
    %% === DoA Estimation Algorithm
    % estimator_angle = DoAEstimator(y_awgn, size(pos_tx,1), lambda, ELEMENT_NUM, element_spacing, sweeping_angle, aoa_true);
    % aoa_rel_est(rx_idx) = estimator_angle.MUSIC().aoa_est;
    % rays_abs{rx_idx} = estimator_coor.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_abs(rx_idx), aoa_rel_est(rx_idx), ABS_ANGLE_LIM);
end


% SSets up a grid for (x,y)
x = linspace(0, 100, 100);
y = linspace(0, 100, 100);
[X, Y] = meshgrid(x, y);

% Calculate theta1 and theta2 from x and y
sin_theta1 = f(X, Y, pos_rx(1,1), pos_rx(1,2), rot_abs(1));
sin_theta2 = f(X, Y, pos_rx(2,1), pos_rx(2,2), rot_abs(2));
% sin_theta1 = f(50, 50, 12.70, 63.24, -57.34)
% sin_theta2 = f(50, 50, 91.34, 9.75, 87.07)

% Calculate the final likelihood function L(x,y)
L = calculateP(sin_theta1, sin_theta2, w{1}, w{2}, ELEMENT_NUM, nPower);

%% === Plotting
fprintf('SNR (dB):\n')
fprintf('%.0f  ', SNR_dB);
fprintf('\n');
fprintf('RX Positions:\n');
for idx = 1:size(pos_rx, 1)
    fprintf('  Rx %d: x = %.2f, y = %.2f\n', idx, pos_rx(idx, 1), pos_rx(idx, 2));
end

fprintf('Absolute RX Rotations (degrees):\n  ');
fprintf('%.2f  ', rot_abs);
fprintf('\n');

fprintf('Absolute RX-TX Angles (degrees):\n  ');
fprintf('%.2f  ', angle_rx_tx_abs);
fprintf('\n');

fprintf('True AoA (degrees):\n  ');
fprintf('%.2f  ', aoa_true);
fprintf('\n');
% 
% fprintf('Estimated AoA (degrees):\n  ');
% fprintf('%.2f  ', aoa_rel_est);
% fprintf('\n');
% Create 3D surface plot
figure;
surf(X, Y, abs(L));
xlabel('x');
ylabel('y');
zlabel('L(x,y)');
title('3D Visualization of L(x,y) = P(\theta_1, \theta_2)');
shading interp;    % Improves plot appearance
colorbar;          % Shows color scale


%--------------------------------------------------------------------------
function sin_theta = f(x, y, x_rx, y_rx, phi)
    % f1 - Computes sin_theta as a function of x and y.
    thetaphi = atan2d(y-y_rx, x-x_rx);
    theta = thetaphi - phi;
    sin_theta = sind(theta);
end

%--------------------------------------------------------------------------
function P = calculateP(sin_theta1, sin_theta2, w1, w2, M, nPower)
    z = [w1.' w2.'].';
    Sigma_z = computeSigmaZ(sin_theta1, sin_theta2, nPower, M);
    P = cell2mat(cellfun(@(Sig) 1/(pi^(2*M)*det(Sig))*exp(-z.'/Sig*z), ...
        Sigma_z, 'UniformOutput', false));
    % (pi^(2*M)*det(Sigma_z))^(-1)*exp(-z.'/Sigma_z*z);
end

function Sigma_z = computeSigmaZ(sin_theta1, sin_theta2, nPower, M)
% Computes the covariance matrix Sigma_z defined as:
%
%   Sigma_z = [ (|gamma1|^2 a1 a1^H + sigma^2 I)   (gamma1 gamma2^* a1 a2^H);
%               (gamma2 gamma1^* a2 a1^H)          (|gamma2|^2 a2 a2^H + sigma^2 I) ]
%   where:
%       - gamma1, gamma2 are the distant coefficients, temporarily set to 1.
%       - a1 = aULA(theta1) and a2 = aULA(theta2) are the steering vectors for a
%         uniform linear array (ULA) of N elements
%   Inputs:
%       theta1 - Angle (in radians) for the first steering vector.
%       theta2 - Angle (in radians) for the second steering vector.
%       nPower - Noise standard deviation (sigma^2 will be used in the matrix).
%       M      - Number of elements in the ULA.
%
%   Output:
%       Sigma_z - The 2N-by-2N covariance matrix.
%
%   Example:
%       theta1 = pi/6;
%       theta2 = pi/4;
%       sigma  = 0.1;
%       N      = 8;
%       Sigma_z = computeSigmaZ(theta1, theta2, sigma, N);

    % Set gamma1 and gamma2 to 1.
    gamma1 = 1;
    gamma2 = 1;
    
    % Compute the steering vectors for the given angles.
    a1 = aULA(sin_theta1, M);
    a2 = aULA(sin_theta2, M);
    
    % Compute the individual blocks.
    block1 = cellfun(@(x) abs(gamma1)^2 * (x * x') + nPower * eye(M), a1, 'UniformOutput', false);
    block2 = cellfun(@(x1, x2) gamma1 * conj(gamma2) * (x1 * x2'), a1, a2, 'UniformOutput', false);
    block3 = cellfun(@(x1, x2) gamma2 * conj(gamma1) * (x1 * x2'), a2, a1, 'UniformOutput', false);
    block4 = cellfun(@(x) abs(gamma2)^2 * (x * x') + nPower * eye(M), a2, 'UniformOutput', false);
    % block1 = abs(gamma1)^2 * (a1 * a1') + nPower^2 * eye(M);
    % block2 = gamma1 * conj(gamma2) * (a1 * a2');
    % block3 = gamma2 * conj(gamma1) * (a2 * a1');
    % block4 = abs(gamma2)^2 * (a2 * a2') + nPower^2 * eye(M);
    
    % Assemble the 2N-by-2N covariance matrix.
    Sigma_z = cellfun(@(b1, b2, b3, b4) [b1, b2; b3, b4], ...
    block1, block2, block3, block4, 'UniformOutput', false);
    % Add a small eye() term to each covariance matrix to avoid singularity
    Sigma_z = cellfun(@(S) S + 1e-9 * eye(size(S)), Sigma_z, 'UniformOutput', false);
end

%--------------------------------------------------------------------------
function a = aULA(sin_theta, M)
% aULA - Generates the steering vector for a ULA of N elements (as a column vector).
%
%   a = aULA(theta, N) returns an N-by-1 steering vector for an array with
%   element spacing corresponding to half the wavelength (hence the factor pi),
%   given by:
%
%       a = [1;
%            exp(-j*pi*sin(theta));
%            exp(-j*pi*2*sin(theta));
%            ...
%            exp(-j*pi*(N-1)*sin(theta)) ]
%
%   Inputs:
%       theta - Angle in radians.
%       N     - Number of array elements.
%
%   Output:
%       a - The N-by-1 steering vector.
    a = arrayfun(@(x) exp(-1i * pi * (0:(M-1))' .* x), ...
    sin_theta, 'UniformOutput', false);
    % [rows, cols] = size(sin_theta);
    % a = zeros(M, rows, cols);
    % for i = 1:M
    %     a(i,:,:) = exp(-1i * pi * (i-1) .* sin_theta);
    % end
end
