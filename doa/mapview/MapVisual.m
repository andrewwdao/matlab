
classdef MapVisual < handle
    properties
        algo_type  % Algorithm type
        tx_pos  % Transmitter position (x, y) in meters
        rx_pos  % Receiver position (x, y) in meters
        area_size  % Area size
        angle_array % Angle array for visualisation  
        powdb_array % Power array in dB for visualisation
        est_aoa % Estimated Angle of Arrival (AoA)
    end
    
    methods
        function obj = MapVisual(algo_type, tx_pos, rx_pos, area_size, angle_array, powdb_array, est_aoa)
            obj.algo_type = algo_type;
            obj.tx_pos = tx_pos;
            obj.rx_pos = rx_pos;
            obj.area_size = area_size;
            obj.angle_array = angle_array;
            obj.powdb_array = powdb_array;
            obj.est_aoa = est_aoa;
        end

        function plotMapOnly(~, area_size, tx_pos, rx_pos, abs_rays, doa_intersect, SHOW_LIMITS)

            % --- Plot the Tx and Rx positions
            plot(tx_pos(:,1), tx_pos(:,2), 'ro', 'MarkerSize', 10, 'LineWidth', 3); hold on;
            text(tx_pos(:,1)+1.5, tx_pos(:,2)-2.5, 'Tx', 'Color', 'red', 'FontSize', 12); hold on;
            plot(rx_pos(:,1), rx_pos(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 3); hold on;
            text(rx_pos(:,1)+1.5, rx_pos(:,2)-2.5, 'Rx', 'Color', 'blue', 'FontSize', 12); hold on;
            
            % --- Plot the rays connecting the Tx and Rx, and the navigation rays
            for i = 1:length(abs_rays)
                fplot(@(x) abs_rays{i}.doa_slope*x + abs_rays{i}.doa_shift, abs_rays{i}.lim, 'r', 'LineWidth', 2); hold on;
            end
            % --- only plotting the DoA intersection if there are 5 input arguments
            if nargin >= 6
                plot(doa_intersect.x, doa_intersect.y, 'ko', 'MarkerSize', 8, 'LineWidth', 2); hold on;
                % --- Plot the navigation rays
                if SHOW_LIMITS
                    for i = 1:length(abs_rays)
                        fplot(@(x) abs_rays{i}.centre_slope*x + abs_rays{i}.centre_shift, abs_rays{i}.lim, 'k--', 'LineWidth', 0.25); hold on;
                        fplot(@(x) abs_rays{i}.cw_slope*x + abs_rays{i}.cw_shift, [abs_rays{i}.lim], 'k--', 'LineWidth', 0.25); hold on;
                        fplot(@(x) abs_rays{i}.ccw_slope*x + abs_rays{i}.ccw_shift, [abs_rays{i}.lim], 'k--', 'LineWidth', 0.25); hold on;
                    end
                end
            
            end

            % Other format settings
            xlabel('X Position (m)'); ylabel('Y Position (m)');
            xlim([0 area_size]); ylim([0 area_size]);
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
            % Plots the normal spectrum of the spatial data.
            % This method plots the spatial spectrum using the angle and power data
            % stored in the object. It creates a subplot in the 2nd row spanning the
            % 3rd and 4th columns, plots the angle vs. power data, adds markers, and
            % labels the axes. It also adds a legend and a title to the plot, and
            % enables the grid.
            %
            % Inputs:
            %    obj - The instance of the class containing the angle and power data.
            %
            % Outputs:
            %    None
            %
            % Example:
            %    obj.plotNormalSpectrum();
            %
            % Other m-files required: None
            % Subfunctions: addMarkers
            % MAT-files required: None
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
            % Plots the spatial spectrum in polar coordinates.
            % This function plots the spatial spectrum data stored in the object in polar coordinates.
            % It normalizes the power spectrum data, plots it using polarplot, and adds markers.
            % The plot is customized with specific angular limits, tick marks, and direction settings.
            %
            % Inputs:
            %    obj - The object containing the spectrum data and angle array.
            %
            % Outputs:
            %    None
            %
            % Example:
            %    obj.plotPolarSpectrum();
            %
            % Other m-files required: None
            % Subfunctions: addPolarMarkers
            % MAT-files required: None
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
        function plotSingle(obj, abs_ray, SHOW_LIMITS)
            figure('Name', 'Map Visualisation with detailed spectrum', 'WindowState', 'maximized'); clf; hold on;
            subplot(2,2,1);
            obj.plotMapOnly(obj.area_size, obj.tx_pos, obj.rx_pos, abs_ray);hold on;
            subplot(2,2,[3,4]);
            obj.plotNormalSpectrum(); hold on;
            subplot(2,2,2);
            obj.plotPolarSpectrum(); hold off;
        end
    end
end