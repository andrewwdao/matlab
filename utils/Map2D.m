
classdef Map2D < handle
    properties
    end
    
    methods
        function obj = Map2D(varargin)
        end

        function plot(obj, pos_tx, pos_rx, rot_rx_abs, area_size, aoa_act, lim_angle, flags, aoa_est_cell)
            num_rx = size(pos_rx, 1);
            num_tx = size(pos_tx, 1);
            rays_abs = cell(num_rx, 1);
            % --- For each Rx
            for rx_idx = 1:num_rx
                % --- Plot the Rx positions
                obj.addMarkers(pos_rx(rx_idx,1), pos_rx(rx_idx,2), ['R',num2str(rx_idx)], 'blue');
                for tx_idx = 1:num_tx
                    % --- Plot the Tx positions
                    obj.addMarkers(pos_tx(tx_idx,1), pos_tx(tx_idx,2), ['T',num2str(tx_idx)], 'k');
                    % Calculate the absolute rays
                    rays_abs{rx_idx} = obj.calAbsRays(pos_rx(rx_idx,:), pos_tx(tx_idx,:), rot_rx_abs(rx_idx), aoa_act(rx_idx), lim_angle);
                    % Draw connecting ray from Tx to Rx
                    plot([pos_tx(tx_idx,1), pos_rx(rx_idx,1)], [pos_tx(tx_idx,2), pos_rx(rx_idx,2)], 'Color', [0 0.5 0], 'LineWidth', 4); hold on;
                end
                % % Determine relative position (from Tx to Rx)
                % delta = pos_rx(rx_idx,:) - pos_tx;
                % if (delta(1) < 0 && delta(2) >= 0) || (delta(1) >= 0 && delta(2) >= 0)
                %     % Quadrant I or II
                %     base_angle = rot_rx_abs(rx_idx);
                % elseif (delta(1) < 0 && delta(2) < 0) || (delta(1) >= 0 && delta(2) < 0)
                %     % Quadrant III or IV
                %     base_angle = rot_rx_abs(rx_idx) + 180;
                % end

                % % Now, sweep the arc from base_angle to base_angle + aoa_act(rx_idx)
                % theta_arc = linspace(base_angle, base_angle + aoa_act(rx_idx), 50); 

                % arc_radius = 5; % Adjust for visual clarity
                % arc_x = pos_rx(rx_idx,1) + arc_radius * cosd(theta_arc);
                % arc_y = pos_rx(rx_idx,2) + arc_radius * sind(theta_arc);
                % plot(arc_x, arc_y, 'm-', 'LineWidth', 2); 
                % hold on;            
            end

            % --- Plot the Estimated AoA intersection rays and point
            if nargin == 9
                rays_est = cell(num_rx, num_tx);
                num_method = size(aoa_est_cell, 1);
                % --- For each method
                for method_idx = 1:num_method
                    % --- For each Tx
                    for tx_idx = 1:num_tx
                        % --- For each Rx
                        for rx_idx = 1:num_rx % Calculate the estimated rays and plot the intersection rays
                           % calculate the correct idx based on the number of Rx and Tx
                            angle_idx = max(rx_idx * (num_rx > num_tx), tx_idx * (num_rx <= num_tx));
                            rays_est{rx_idx, tx_idx} = obj.calAbsRays(pos_rx(rx_idx,:), pos_tx(tx_idx,:), rot_rx_abs(rx_idx), aoa_est_cell{method_idx, angle_idx});
                            fplot(@(x) rays_est{rx_idx, tx_idx}.doa_slope*x + rays_est{rx_idx, tx_idx}.doa_shift, rays_est{rx_idx, tx_idx}.lim, 'r', 'LineWidth', 2); hold on;
                        end
                        
                        % Only plot the intersection point if there are more than 1 Rx
                        if num_rx > 1
                            for rx_idx = 1:2:num_rx 
                                % Calculate the intersection point of the two rays
                                doa_intersect = obj.calDoAIntersect(rays_est{rx_idx, tx_idx}, rays_est{rx_idx+1, tx_idx});
                                % Plot the intersection point
                                plot(doa_intersect.x, doa_intersect.y, 'mo', 'MarkerSize', 8, 'LineWidth', 3); hold on;
                            end
                        end
                    end
                end
            end

            show_limits = flags(1);
            show_extra = flags(2);
            % --- Plot the navigation rays
            if show_limits
                for rx_idx = 1:num_rx
                    fplot(@(x) rays_abs{rx_idx}.cw_slope*x + rays_abs{rx_idx}.cw_shift, [rays_abs{rx_idx}.lim], 'k--', 'LineWidth', 0.25); hold on;
                    fplot(@(x) rays_abs{rx_idx}.centre_slope*x + rays_abs{rx_idx}.centre_shift, rays_abs{rx_idx}.lim, 'k--', 'LineWidth', 0.25); hold on;
                    fplot(@(x) rays_abs{rx_idx}.ccw_slope*x + rays_abs{rx_idx}.ccw_shift, [rays_abs{rx_idx}.lim], 'k--', 'LineWidth', 0.25); hold on;
                end
            end

            if show_extra
                % --- Plot the Tx positions
                for tx_idx = 1:num_tx
                    obj.addMarkers(pos_tx(tx_idx,1), pos_tx(tx_idx,2), ['T',num2str(tx_idx)], 'k');
                end
                % --- Plot the Rx positions
                for rx_idx = 1:num_rx
                    % add the true aoa in the plot
                    text(pos_rx(rx_idx,1)+1.5, pos_rx(rx_idx,2)-7.5, sprintf('AoA: %.2f°', aoa_act(rx_idx)), 'Color', 'blue', 'FontSize', 12); hold on;
                    % add the estimated aoa in the plot
                    text(pos_rx(rx_idx,1)+1.5, pos_rx(rx_idx,2)-12.5, sprintf('Est. AoA: %.2f°', aoa_est_cell{1, rx_idx}), 'Color', 'blue', 'FontSize', 12); hold on;
                
                end
            end

            % Other format settings
            xlabel('X Position (m)'); ylabel('Y Position (m)');
            xlim([0 area_size]); ylim([0 area_size]);
            title('Map Visualisation'); grid on; hold off;
        end

        function plotDetailed(obj, pos_tx, pos_rx, rot_rx_abs, area_size, aoa_act, lim_angle, flags, angle_array, powdb_cell, method_list, aoa_est_cell)
            figure('Name', 'Map and Spectrum Visualisation', 'WindowState', 'maximized'); clf; hold on;
            subplot(2,2,1);
            obj.plot(pos_tx, pos_rx, rot_rx_abs, area_size, aoa_act, lim_angle, flags, aoa_est_cell);hold on;
            subplot(2,2,[3,4]);
            obj.plotSpectrum(angle_array, powdb_cell, method_list); hold on;
            subplot(2,2,2);
            obj.plotPolarSpectrum(angle_array, powdb_cell, method_list); hold off;
        end

        function plotSpectrum(obj, angle_array, powdb_cell, method_list, aoa_est_cell)
            method_num = size(powdb_cell, 1);
            tx_num = size(powdb_cell, 2);
            pl = gobjects(1, method_num);
            for method_idx = 1:method_num % only need to use the first spectrum as they would be the same for all
                combinedSpectrum = zeros(size(angle_array));
                % Get the combined spectrum first
                for tx_idx = 1:tx_num
                    combinedSpectrum = combinedSpectrum + powdb_cell{method_idx,tx_idx};
                end
                pl(method_idx) = plot(angle_array, combinedSpectrum, 'LineWidth', 2); hold on;
                % Add markers for the estimated AoA if required
                if nargin == 5
                    for tx_idx = 1: tx_num
                        [est_pow, ~] = maxk(combinedSpectrum, tx_num); % Get the peaks to visualise the AoA on the spectrum plot 
                        obj.addMarkers(aoa_est_cell{method_idx, tx_idx}, est_pow(tx_idx), ['DoA:', num2str(aoa_est_cell{method_idx, tx_idx}), '°; P:', num2str(est_pow(tx_idx)), 'dB'], 'red');
                    end
                end
            end
            xlabel('Angle (degrees)');
            ylabel('Power (dB)');
            legend(pl, method_list, 'AutoUpdate', 'off');
            title('Spatial Spectrum'); grid on; hold off;
        end

        function plotPolarSpectrum(obj, angle_array, powdb_cell, method_list, aoa_est_cell)
            method_num = size(powdb_cell, 1);
            tx_num = size(powdb_cell, 2);
            pl = gobjects(1, method_num);
            for method_idx = 1:method_num % only need to use the first spectrum as they would be the same for all
                powdb_cell_normalised = cellfun(@(x) x - min(x), powdb_cell, 'UniformOutput', false); % Normalize the spectrum data
                combinedSpectrum = zeros(size(angle_array));
                for tx_idx = 1:tx_num
                    combinedSpectrum = combinedSpectrum + powdb_cell_normalised{method_idx,tx_idx};
                end
                pl(method_idx) = polarplot(deg2rad(angle_array), combinedSpectrum, '-', 'LineWidth', 2); hold on;
                % Add markers for the estimated AoA if required
                if nargin == 5
                    for tx_idx = 1: tx_num
                        [est_pow, ~] = maxk(combinedSpectrum, tx_num); % Get the peaks to visualise the AoA on the spectrum plot 
                        obj.addPolarMarkers(aoa_est_cell{method_idx, tx_idx}, est_pow(tx_idx), ['DoA:', num2str(aoa_est_cell{method_idx, tx_idx}), '°; P:', num2str(est_pow(tx_idx)), 'dB'], 'red');
                    end
                end
            end
            ax = gca;
            ax.RTickLabel = '';
            ax.ThetaLim = [min(angle_array) max(angle_array)];
            ax.ThetaTick = min(angle_array):15:max(angle_array);
            ax.ThetaZeroLocation = 'right';  % 0 degrees at the right
            ax.ThetaDir = 'counterclockwise';  % Counterclockwise direction
            legend(pl, method_list, 'AutoUpdate', 'off');
            title('Spatial Spectrum (Polar)');
        end

        function addMarkers(~, x, y, text_str, color)
            %   Plots a marker at the given (x, y) coordinates with an associated text label.
            %
            % Inputs:
            %   x         - The x-coordinate for the marker.
            %   y         - The y-coordinate for the marker.
            %   text_str  - The text string to display near the marker.
            %   color     - The color specification for both the marker and the text.
            %
            % Notes:
            %   - A circle marker is plotted at (x, y) with a specified marker size and line width.
            %   - A text label is added with an offset from the (x, y) coordinates.
            %   - 'hold on' is called to retain the current plot when adding new elements.
            plot(x, y, 'Marker', 'o', 'Color', color, 'MarkerSize', 10, 'LineWidth', 3, 'LineStyle', 'none'); hold on;
            text(x+1.5, y-2.5, text_str, 'Color', color, 'FontSize', 12); hold on;
        end

        function addPolarMarkers(~, angle, radius, text_str, color)
            polarplot(deg2rad(angle), radius, 'Marker', 'o', 'Color', color, 'MarkerSize', 10, 'LineWidth', 3, 'LineStyle', 'none'); hold on;
            text(deg2rad(angle), radius*1.3, text_str, 'Color', color, 'FontSize', 12); hold on;
        end

        function abs_ray = calAbsRays(obj, pos_abs_rx, pos_abs_tx, rot_abs_rx, aoa_act, rot_lim)
            % Calculate the absolute ray parameters (slope and doa_shift)
            % for a line connecting the Tx and Rx in the 2D plane.
            %
            % This is an absolute slope and doa_shift added to
            % a relative angle of arrival for an RX at the origin (0,0) with no doa_shift.
            % 
            % The function also determines the limits of the ray based on the quadrant
            % in which the Rx and Tx are located.
            %
            % Inputs:
            %    pos_abs_rx - A 1x2 vector representing the absolute position of the Rx in the 2D plane [x, y].
            %    pos_abs_tx - A 1x2 vector representing the absolute position of the Tx in the 2D plane [x, y].
            %    rot_abs_rx - A scalar representing the absolute rotation of the Rx in degrees.
            %    aoa_act - A scalar representing the true relative angle of arrival from TX to RX in degrees.
            %    rot_lim (optional) - A scalar representing the limit of the rotation in degrees.
            %
            % Outputs:
            %    abs_ray    - A structure containing the following fields:
            %                 centre_slope - The slope of the line of the absolute rotation of the Rx. to the x-axis.
            %                 cw_slope     - The slope of the line of the clock-wise absolute rotation of the Rx to the x-axis.
            %                 ccw_slope    - The slope of the line of the counter-clock-wise absolute rotation of the Rx to the x-axis.
            %                 doa_slope    - The slope of the line of the absolute angle of arrival of the Rx to the Tx.
            %                 doa_shift - The absolute doa_shift of the lines that go through RX.
            %                 lim   - A 1x2 vector representing the limits of the ray in the 2D plane that is used for plotting.
            %
            % Example:
            %    pos_abs_rx = [100, 200];
            %    pos_abs_tx = [300, 400];
            %    rot_abs_rx = 45;
            %    aoa_act = 30;
            abs_ray.doa_slope = tand(rot_abs_rx + aoa_act);
            abs_ray.doa_shift = pos_abs_rx(2) - abs_ray.doa_slope * pos_abs_rx(1);
            abs_ray.lim = obj.calRayPlotLim(pos_abs_rx, pos_abs_tx, abs_ray.doa_slope);
            if nargin == 6
                abs_ray.centre_slope = tand(rot_abs_rx);
                abs_ray.centre_shift = pos_abs_rx(2) - abs_ray.centre_slope * pos_abs_rx(1);
                abs_ray.cw_slope = tand(rot_abs_rx - rot_lim);
                abs_ray.cw_shift = pos_abs_rx(2) - abs_ray.cw_slope * pos_abs_rx(1);
                abs_ray.ccw_slope = tand(rot_abs_rx + rot_lim);
                abs_ray.ccw_shift = pos_abs_rx(2) - abs_ray.ccw_slope * pos_abs_rx(1);
            end
        end

        function ray_plot_lim = calRayPlotLim(~, pos_abs_rx, pos_abs_tx, slope)
            % Calculate the limits of the ray plot based on the quadrant in which the Rx and Tx are located.
            % 
            % Inputs:
            %    pos_abs_rx - A 1x2 vector representing the absolute position of the Rx in the 2D plane [x, y].
            %    pos_abs_tx - A 1x2 vector representing the absolute position of the Tx in the 2D plane [x, y].
            %    slope      - A scalar representing the slope of the line connecting the Tx and Rx.
            % 
            % Outputs:
            %    ray_plot_lim - A 1x2 vector representing the limits of the ray plot in the 2D plane.
            % 
            % Example:
            %    pos_abs_rx = [100, 200];
            %    pos_abs_tx = [300, 400];
            %    slope = 0.5;
            if (pos_abs_rx(1) < pos_abs_tx(1) && pos_abs_rx(2) < pos_abs_tx(2)) || (pos_abs_rx(1) > pos_abs_tx(1) && pos_abs_rx(2) < pos_abs_tx(2))
                lim_range = [0, 10000]; % Third quadrant or first quadrant
            else
                lim_range = [-10000, 0]; % Second quadrant or fourth quadrant
            end
            ray_plot_lim = sort(slope * lim_range + pos_abs_rx(1), "ascend"); % Absolute ray limit from a relative x limit from 0 to 100
        end

        function result = calDoAIntersect(~, abs_ray1, abs_ray2)
            % Calculate the intersection point of two lines in the 2D plane.
            % 
            % Inputs:
            %    abs_ray1 - The first line in the 2D plane.
            %    abs_ray2 - The second line in the 2D plane.
            %
            % Outputs:
            %    result - The intersection point of the two lines.
            %
            % Other m-files required: None
            % Subfunctions: None
            % MAT-files required: None
            result.x = (abs_ray2.doa_shift - abs_ray1.doa_shift) / (abs_ray1.doa_slope - abs_ray2.doa_slope);
            result.y = abs_ray1.doa_slope * result.x + abs_ray1.doa_shift;
        end
    end
end