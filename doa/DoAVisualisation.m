
classdef DoAVisualisation < handle
    properties
        algo_type
        tx_pos
        rx_pos
        area_size
        angle_array
        powdb_array
        est_aoa
    end
    
    methods
        function obj = DoAVisualisation(algo_type, tx_pos, rx_pos, area_size, angle_array, powdb_array, est_aoa)
            obj.algo_type = algo_type;
            obj.tx_pos = tx_pos;
            obj.rx_pos = rx_pos;
            obj.area_size = area_size;
            obj.angle_array = angle_array;
            obj.powdb_array = powdb_array;
            obj.est_aoa = est_aoa;
        end
        function plotMap(obj)
            figure('Name', 'Spatial Spectrum', 'WindowState', 'maximized'); clf;
            subplot(2,2,1); hold on;
            plot(obj.tx_pos(:,1), obj.tx_pos(:,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            text(obj.tx_pos(:,1) + 2, obj.tx_pos(:,2), 'Tx', 'Color', 'red', 'FontSize', 12);
            plot(obj.rx_pos(:,1), obj.rx_pos(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
            text(obj.rx_pos(:,1) + 2, obj.rx_pos(:,2), 'Rx', 'Color', 'blue', 'FontSize', 12);
            xlim([0 obj.area_size]); ylim([0 obj.area_size]);
            xlabel('X Position (m)'); ylabel('Y Position (m)');
            title('Map Visualisation'); grid on; hold off;
        end
        function addMarkers(obj)
            [est_pow, ~] = maxk(obj.powdb_array, length(obj.est_aoa)); % Get the selected maximum pow_array_db 
            for i = 1:length(obj.est_aoa)
                plot(obj.est_aoa(i), est_pow(i), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2);
                text(obj.est_aoa(i)-5, (est_pow(i)+0.5).*1.15, ['DoA:', num2str(obj.est_aoa(i)), '°; P:', num2str(est_pow(i)), 'dB'], 'LineWidth', 2);
            end
        end
        function plotNormalSpectrum(obj)
            subplot(2,2,[3,4]);
            plot(obj.angle_array, obj.powdb_array, 'LineWidth', 2); hold on;
            obj.addMarkers();
            xlabel('Angle (degrees)');
            ylabel('Power (dB)');
            legend(obj.algo_type, 'AutoUpdate', 'off');
            title('Spatial Spectrum'); grid on; hold off;
        end
        function addPolarMarkers(obj, powdb_array_normalized)
            [est_pow, est_idx] = maxk(obj.powdb_array, length(obj.est_aoa)); % Get the selected maximum
            for i = 1:length(obj.est_aoa)
                polarplot(deg2rad(obj.est_aoa(i)), powdb_array_normalized(est_idx(i)), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2);
                text(deg2rad(obj.est_aoa(i)), powdb_array_normalized(est_idx(i))*1.15, ['DoA:', num2str(obj.est_aoa(i)), '°; P:', num2str(est_pow(i)), 'dB'], 'LineWidth', 2);
            end
        end
        function plotPolarSpectrum(obj)
            subplot(2,2,2);
            powdb_array_normalized = obj.powdb_array - min(obj.powdb_array); % Normalize the spectrum data
            polarplot(deg2rad(obj.angle_array), powdb_array_normalized, '-', 'LineWidth', 2); hold on;
            obj.addPolarMarkers(powdb_array_normalized);
            hold off;
            ax = gca;
            ax.RTickLabel = '';
            ax.ThetaLim = [min(obj.angle_array) max(obj.angle_array)];
            ax.ThetaTick = min(obj.angle_array):15:max(obj.angle_array);
            ax.ThetaZeroLocation = 'right';  % 0 degrees at the right
            ax.ThetaDir = 'counterclockwise';  % Counterclockwise direction
            title('Spatial Spectrum (Polar)');
        end
        function plot(obj)
            obj.plotMap();hold on;
            obj.plotNormalSpectrum(); hold on;
            obj.plotPolarSpectrum(); hold off;
        end
    end
end