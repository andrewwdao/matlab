%Initialisation and visualise Tx/Rx on a map
clear; clc; close all;
% Constants
c = physconst('LightSpeed');
fc = 2.4e9; % Operating frequency in Hz - 2.4 GHz
lambda = c / fc;
% Area size and positions
area_size = 100;         % 100x100 meter area
tx_position = [80, 80];  % Transmitter position (x, y) in meters
rx_position = [20, 20];  % Receiver position (x, y) in meters
% Antenna Array Configuration
Array = phased.ULA('NumElements', 64, 'ElementSpacing', 0.5);  % 8x8 antenna array
pos = getElementPosition(Array)/lambda;  % Element positions
% Signal and noise parameters
Nsamp = 1000;
nPower = 0.1; % W
% True Angle of Arrival (AoA)
true_aoa = atan2d(tx_position(2) - rx_position(2), tx_position(1) - rx_position(1));

% Generate signal received at Rx
rs=rng(2007); % initialize the random number generator in MATLAB to a specific seed value
signal = sensorsig(pos, Nsamp, true_aoa, nPower);

% MVDR Estimator
mvdrspatialspect = phased.MVDREstimator('SensorArray',Array,...
    'OperatingFrequency',fc,'ScanAngles',-90:90,...
    'DOAOutputPort',true,'NumSignals',2);
% Estimate DoA using MVDR
[~, estimated_aoa] = mvdrspatialspect(signal);

% Get the spectrum data
% scan_angles = mvdr_estimator.AzimuthScanAngles;  % Retrieve the scan angles
scan_angles = mvdrspatialspect.ScanAngles;  % Retrieve the scan angles
ymvdr = mvdrspatialspect(signal);           % Get the spectrum data
% Convert spectrum data to dB
ymvdr_dB = 20 * log10(ymvdr);
% Normalize the spectrum data (optional: you can skip this if you want absolute values)
ymvdr_dB_normalized = ymvdr_dB - min(ymvdr_dB);  % Normalize to the peak


fprintf('\nTrue Angles of Arrival: %f', true_aoa)
fprintf('\nEstimated Angle of Arrival: %f', estimated_aoa)
% --- Plotting
% -- Map
figure;
plot(tx_position(1), tx_position(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
plot(rx_position(1), rx_position(2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
text(tx_position(1) + 2, tx_position(2), 'Tx', 'Color', 'red', 'FontSize', 12);
text(rx_position(1) + 2, rx_position(2), 'Rx', 'Color', 'blue', 'FontSize', 12);
xlim([0 area_size]);
ylim([0 area_size]);
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Map with Tx and Rx Positions');
legend('Tx Position', 'Rx Position');
grid on;
hold off;
% -- Normal spectrum
figure; title('DoA Estimation using MVDR with ULA');
plotSpectrum(mvdrspatialspect);
title('MVDR Spatial Spectrum');
% -- Normalized spectrum on a polar plot
figure; clf;
% Create a new figure for the polar plot or set it to a new axes
% axes('Position', [0.35 0.35 0.3 0.3]); % Positioning the radar plot on the map
polarplot(deg2rad(scan_angles), ymvdr_dB_normalized, 'r-', 'LineWidth', 2);
% --- Customize the polar plot
ax = gca;
ax.RTickLabel = '';
ax.ThetaLim = [0 360];
ax.ThetaTick = 0:15:360;
ax.ThetaZeroLocation = 'top';  % 0 degrees at the right
ax.ThetaDir = 'clockwise';  % Counterclockwise direction
% --- Mark the estimated DoA
hold on;
% Angular resolution for scanning
theta = -90:1:90; % Scanning angles in degrees
[max_pow, max_idx] = max(ymvdr_dB);
polarplot(deg2rad(estimated_aoa(1)), ymvdr_dB_normalized(max_idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(deg2rad(estimated_aoa(1))-pi/12, ymvdr_dB_normalized(max_idx)-3, ...
    ['DoA:', num2str(estimated_aoa(1)), '°; P:', num2str(max_pow), 'dB'], 'Color', 'red');
hold off;
title('MVDR Spatial Spectrum (Polar)');

