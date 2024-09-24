%Initialisation and visualise Tx/Rx on a map
clear; clc; close all;

% Area size and positions
area_size = 100;    % 100x100 meter area
tx_pos = [60,80;];  % Transmitter position (x, y) in meters
rx_pos = [10,50;];  % Receiver position (x, y) in meters
element_num = 4;  % Number of elements in the ULA
% Signal and noise parameters
Nsamp = 1000;
nPower_db = 10; % White noise power (dB)
% True Angle of Arrival (AoA)
true_aoa = zeros(size(rx_pos, 1), size(tx_pos, 1));
for i = 1:size(rx_pos, 1)
    for j = 1:size(tx_pos, 1)
        true_aoa(i,j) = atan2d(tx_pos(j,2) - rx_pos(i,2), tx_pos(j,1) - rx_pos(i,1));
    end
end
disp('---------------------------------- True Angles of Arrival:');
disp(array2table(...
    true_aoa, ...% table data
    'RowNames', cellstr(strcat('RX', num2str((1:size(rx_pos, 1))'))), ...
    'VariableNames', cellstr(strcat('TX', num2str((1:size(tx_pos, 1))')))));



% ====================== Generate signal received at Rx =======================
% Area size and positions
% area_size = 100;    % 100x100 meter area
% tx_pos = [60,80;];  % Transmitter position (x, y) in meters
% rx_pos = [10,50;];  % Receiver position (x, y) in meters
% element_num = 4;  % Number of elements in the ULA
% Constants
c = 299792458; % physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength
t = (0:1e-6:1e-3);  % Time vector for the signal
N = element_num;  % Number of antenna elements in the ULA
element_spacing = 0.5 * lambda;  % Element spacing (ULA)
P_t = 1;                       % Transmit signal power
% --- Calculate the Angle of Arrival (AoA)
theta = atan2d(tx_pos(2) - rx_pos(2), tx_pos(1) - rx_pos(1));  % AoA in degrees
rs=rng(2007); % initialize the random number generator to a specific seed value
% --- Generate transmitted signal
% Transmitted signal (simple sinusoidal signal)
s_t = sqrt(P_t) * exp(1j * 2 * pi * fc * t);  % Complex sinusoid
% --- Transmitted through a Friis free-space model
% Calculate the Euclidean distance
distance = sqrt((tx_pos(1) - rx_pos(1))^2 + (tx_pos(2) - rx_pos(2))^2);
channel = ChannelModel(lambda, distance, 1);
y_t = channel.FriisModel(s_t);  % Received signal at the receiver

% --- MUSIC Algorithm
% --- Calculate the covariance matrix
R = y_t * y_t' / size(y_t, 2);
% Perform eigenvalue decomposition
[eigenvectors, eigenvalues] = eig(R);
% Sort eigenvalues and eigenvectors
[eigenvalues, idx] = sort(diag(eigenvalues), 'descend');
eigenvectors = eigenvectors(:, idx);
% Determine the noise subspace
num_signals = 1; % Number of signals (assuming 1 for simplicity)
noise_subspace = eigenvectors(:, num_signals+1:end);
% Compute the MUSIC spectrum
angles = -90:1:90; % Angle range for MUSIC spectrum
music_spectrum = zeros(size(angles));

for i = 1:length(angles)
    steering_vector = exp(-1j * 2 * pi * element_spacing * (0:N-1)' * sind(angles(i)) / lambda);
    music_spectrum(i) = 1 / (steering_vector' * (noise_subspace * noise_subspace') * steering_vector +eps(1)); % add a small positive constant to prevent division by zero. 9.44 in [1]
end

% Convert MUSIC spectrum to dB scale
music_spectrum_dB = 10 * log10(abs(music_spectrum));

% Plot the MUSIC spectrum
figure;
plot(angles, music_spectrum_dB);
title('MUSIC Spectrum');
xlabel('Angle (degrees)');
ylabel('Spectrum (dB)');
grid on;

% Find the peaks in the MUSIC spectrum
[~, peak_indices] = findpeaks(music_spectrum_dB);
estimated_aoa = angles(peak_indices);

% Display the estimated angles of arrival
disp('Estimated Angles of Arrival (AoA):');
disp(peak_indices);
% --- Steering vector for the incoming signal
steering_vec = exp(-1j * 2 * pi * distance * (0:N-1)' * sind(theta) / lambda);
received_signal = y_t' * steering_vec';  % Received signal at the receiver





% Define the Antenna Array Configuration
Array = phased.ULA('NumElements', element_num, 'ElementSpacing', 0.5*lambda);  % wavelength spacing
% POS = getElementPosition(H) returns the element positions of
%   the ULA H. POS is a 3xN matrix where N is the number of
%   elements in H. Each column of POS defines the position, in the
%   form of [x; y; z] (in meters), of an element in the local
%   coordinate system.
pos = getElementPosition(Array)/lambda;  % Element positions
% --------- Simulate received signal at the sensor array
%   [X,RT,R] = sensorsig(POS,NS,ANG,NCOV,SCOV)
%   - POS represents the locations of elements in the sensor array, specified
%   in the unit of signal wavelength. All elements in the sensor array are
%   assumed to be isotropic.
%   - NS is the number of snapshots.
%   - ANG is the directions of the incoming signals in degrees.
%   If ANG is a 2xM matrix, each column specifies the direction in
%   the space in [azimuth; elevation] form (in degrees)
%   - NCOV is the noise noise power (in watts).
%   If NCOV is a nonnegative scalar, it representsthe white noise across all sensor elements.
%   If NCOV is a 1xN vector, it represents the noise power (in watts) of each receiving sensor.
%   If NCOV is an NxN matrix, then it represents the covariance matrix for the noise across all sensor elements.
%   - SCOV is the signal power (in watts).
%   If SCOV is a scalar, it represents the power (in watts) of the incoming signals.
%   All incoming signals are assumed to be uncorrelated and share the same power level.
%   If SCOV is a 1xM matrix, then it represents the power of individual incoming
%   signals. However, all incoming signals are still uncorrelated to each other.
%   If SCOV is an MxM matrix, then it represents the covariance matrix for all incoming signals.
%   - X is the received signal in an NSxN matrix. Each column of X represents
%   the received signal at the corresponding element.
%   - RT returns the theoretical covariance matrix of the received signal in an NxN matrix.
%   - R returns the sample covariance matrix of the received signal in an NxN matrix.
%   R is derived from X.
%   The input signals are assumed to be constant-modulus signals with random phases.
[signal, ~, R] = sensorsig(pos, Nsamp, true_aoa, db2pow(nPower_db));
noise = 0.1*(randn(size(signal))+1i*randn(size(signal)));
% --- Conventional Beamforming
estimator = phased.MUSICEstimator(...
    'SensorArray', Array,...
    'PropagationSpeed', c, 'OperatingFrequency', fc, 'ScanAngles', -90:0.5:90,...
    'DOAOutputPort', true, 'NumSignalsSource', 'Property', 'NumSignals', size(tx_pos, 1));
% Estimate DoA using Conventional
[ypow, est_aoa] = estimator(signal+noise); % Get the spectrum data and the estimated AoA
ypow_dB = 20*log10(ypow) - max(20*log10(ypow)); % Convert spectrum data to dB
[max_pow, ~] = maxk(ypow_dB, length(est_aoa)); % Get the selected maximum power and its index
disp('---------------------------------- Estimated Barlett AoA:');
disp(array2table(...
    est_aoa, ...% table data
    'RowNames', cellstr(strcat('RX', num2str((1:size(rx_pos, 1))'))), ...
    'VariableNames', cellstr(strcat('TX', num2str((1:size(tx_pos, 1))')))));

% % --- Plotting
vis = DoAVisualisation(...
    "MUSIC", tx_pos, rx_pos, area_size, ...
    estimator.ScanAngles, ypow_dB, ...
    est_aoa);
vis.plot();
