%Initialisation and visualise Tx/Rx on a map
clear; clc; close all;

% Area size and positions
area_size = 100;     % 100x100 meter area
tx_pos = [80,80; 75,63; 82,96;];  % Transmitter position (x, y) in meters
rx_pos = [20, 20];  % Receiver position (x, y) in meters
element_num = 64;     % Number of elements in the ULA
% Signal and noise parameters
Nsamp = 1000;
nPower_db = 10; % White noise power (dB)
% Constants
c = physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz)
lambda = c / fc; % Wavelength

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
% Antenna Array Configuration
Array = phased.ULA('NumElements', element_num, 'ElementSpacing', 0.5*lambda);  % antenna array with 0.5 wavelength spacing


% Generate signal received at Rx
rs=rng(2007); % initialize the random number generator in MATLAB to a specific seed value
pos = getElementPosition(Array)/lambda;  % Element positions
signal = sensorsig(pos, Nsamp, true_aoa, db2pow(nPower_db)); % Simulate received signal at sensor array

% --- Conventional Beamforming
bartlettspect = phased.BeamscanEstimator(...
    'SensorArray',Array,...
    'PropagationSpeed', c, 'OperatingFrequency',fc,'ScanAngles',-90:90,...
    'DOAOutputPort', true, 'NumSignals', size(tx_pos, 1));
% Estimate DoA using Conventional
[yconv, est_aoa_conv] = bartlettspect(signal); % Get the spectrum data and the estimated AoA
yconv_dB = 20*log10(yconv) - max(20*log10(yconv)); % Convert spectrum data to dB
[max_pow_conv, max_idx_conv] = maxk(yconv_dB, length(est_aoa_conv)); % Get the selected maximum power and its index
disp('---------------------------------- Estimated Barlett AoA:');
disp(array2table(...
    est_aoa_conv, ...% table data
    'RowNames', cellstr(strcat('RX', num2str((1:size(rx_pos, 1))'))), ...
    'VariableNames', cellstr(strcat('TX', num2str((1:size(tx_pos, 1))')))));
% --- MVDR Estimator - Minimum Variance Distortion-less Response
% The MVDR algorithm's improved resolution comes with a price.
% The MVDR is more sensitive to sensor position errors.
% In circumstances where sensor positions are inaccurate, MVDR could produce a worse spatial spectrum than beamscan.
% Moreover, if we further reduce the difference of two signal directions to a level that is smaller than the beamwidth of an MVDR beam, the MVDR estimator will also fail.
mvdrspatialspect = phased.MVDREstimator(...
    'SensorArray', Array,...
    'PropagationSpeed', c, 'OperatingFrequency', fc, 'ScanAngles', -90:90,...
    'DOAOutputPort', true, 'NumSignals', size(tx_pos, 1));
% Estimate DoA using MVDR
[ymvdr, est_aoa_mvdr] = mvdrspatialspect(signal); % Get the spectrum data and the estimated AoA
ymvdr_dB = 20*log10(ymvdr) - max(20*log10(ymvdr)); % Convert spectrum data to dB
[max_pow_mvdr, max_idx_mvdr] = maxk(ymvdr_dB, length(est_aoa_mvdr)); % Get the selected maximum power and its index
disp('---------------------------------- Estimated MVDR AoA:');
disp(array2table(...
    est_aoa_mvdr, ...% table data
    'RowNames', cellstr(strcat('RX', num2str((1:size(rx_pos, 1))'))), ...
    'VariableNames', cellstr(strcat('TX', num2str((1:size(tx_pos, 1))')))));
% --- MUSIC Estimator - Multiple Signal Classification
% MUSIC provides better spatial resolution than MVDR.
% However, like MVDR, is sensitive to sensor position errors.
% In addition, the number of sources must be known or accurately estimated.
% When the number of sources specified is incorrect, MVDR and Beamscan may simply return insignificant peaks from the correct spatial spectrum.
% In contrast, the MUSIC spatial spectrum itself may be inaccurate when the number of sources is not specified correctly.
% In addition, the amplitudes of MUSIC spectral peaks cannot be interpreted as the power of the sources.
musicspatialspect = phased.MUSICEstimator(...
    'SensorArray', Array,...
    'PropagationSpeed', c, 'OperatingFrequency', fc, 'ScanAngles', -90:90,...
    'DOAOutputPort', true, 'NumSignalsSource', 'Property', 'NumSignals', size(tx_pos, 1));
[ymusic, est_aoa_music] = musicspatialspect(signal); % Get the spectrum data and the estimated AoA
ymusic_dB = 20*log10(ymusic) - max(20*log10(ymusic)); % Convert spectrum data to dB
[max_pow_music, max_idx_music] = maxk(ymusic_dB, length(est_aoa_music)); % Get the selected maximum power and its index
disp('---------------------------------- Estimated MUSIC AoA:');
disp(array2table(... 
    est_aoa_music, ...% table data
    'RowNames', cellstr(strcat('RX', num2str((1:size(rx_pos, 1))'))), ...
    'VariableNames', cellstr(strcat('TX', num2str((1:size(tx_pos, 1))')))));

% --- Plotting
figure('Name', 'Spatial Spectrum', 'WindowState', 'maximized'); clf;
% -- Map
subplot(2,2,1); hold on;
plot(tx_pos(:,1), tx_pos(:,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(tx_pos(:,1) + 2, tx_pos(:,2), 'Tx', 'Color', 'red', 'FontSize', 12);
plot(rx_pos(:,1), rx_pos(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
text(rx_pos(:,1) + 2, rx_pos(:,2), 'Rx', 'Color', 'blue', 'FontSize', 12);
xlim([0 area_size]);
ylim([0 area_size]);
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Map with Tx and Rx Positions');
legend('Tx Position', 'Rx Position');
grid on;
hold off;
% axes('Position', [0.35 0.35 0.3 0.3]); % Positioning the radar plot on the map
% -- Normal spectrum
subplot(2,2,[3,4]);
% plotSpectrum(mvdrspatialspect);
% Plot spectra in dB normalized to 0 dB at the peak
plot(...
    bartlettspect.ScanAngles, yconv_dB, ...
    mvdrspatialspect.ScanAngles, ymvdr_dB, ...
    musicspatialspect.ScanAngles, ymusic_dB, ...
    'LineWidth', 2);
xlabel('Angle (degrees)');  % Default for ULA - xlabel('Elevation Angle (degrees)'); for URA
ylabel('Power (dB)');
legend('Conventional', 'MVDR', 'MUSIC', 'AutoUpdate', 'off');
grid on;
title('Spatial Spectrum');
% Add markers for the top spectrum peaks
hold on;
for i = 1:length(est_aoa_mvdr)
    plot(mvdrspatialspect.ScanAngles(max_idx_mvdr(i)), ymvdr_dB(max_idx_mvdr(i)), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
    text(mvdrspatialspect.ScanAngles(max_idx_mvdr(i))-5, ymvdr_dB(max_idx_mvdr(i))+1, ...
        ['DoA:', num2str(est_aoa_mvdr(i)), '°; P:', num2str(max_pow_mvdr(i)), 'dB'], 'Color', 'blue', 'LineWidth', 2);
end

for i = 1:length(est_aoa_music)
    plot(musicspatialspect.ScanAngles(max_idx_music(i)), ymusic_dB(max_idx_music(i)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(musicspatialspect.ScanAngles(max_idx_music(i))-5, ymusic_dB(max_idx_music(i))+1, ...
        ['DoA:', num2str(est_aoa_music(i)), '°; P:', num2str(max_pow_music(i)), 'dB'], 'Color', 'red', 'LineWidth', 2);
end
hold off;
% -- Normalized spectrum on a polar plot
subplot(2,2,2);
% scan_angles = mvdr_estimator.AzimuthScanAngles;  % Retrieve the scan angles
ymvdr_dB_normalized = ymvdr_dB - min(ymvdr_dB);  % Normalize the spectrum data to 0 dB at the peak
% compress the spectrum data to fit the polar plot with a marker at the estimated DoA
[max_pow, max_idx] = maxk(ymvdr_dB, length(est_aoa_mvdr));
polarplot(deg2rad(mvdrspatialspect.ScanAngles), ymvdr_dB_normalized, 'r-', 'LineWidth', 2); hold on;
for i = 1:length(est_aoa_mvdr)
    polarplot(deg2rad(est_aoa_mvdr(i)), ymvdr_dB_normalized(max_idx(i)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(deg2rad(est_aoa_mvdr(i)), ymvdr_dB_normalized(max_idx(i))+2, ...
        ['DoA:', num2str(est_aoa_mvdr(i)), '°; P:', num2str(max_pow(i)), 'dB'], 'Color', 'red', 'LineWidth', 2);
end
hold off;
% --- Customize the polar plot
ax = gca;
ax.RTickLabel = '';
ax.ThetaLim = [0 360];
ax.ThetaTick = 0:15:360;
ax.ThetaZeroLocation = 'right';  % 0 degrees at the right
ax.ThetaDir = 'counterclockwise';  % Counterclockwise direction
title('Spatial Spectrum (Polar)');
