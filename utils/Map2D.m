
classdef Map2D < handle
    properties
    end
    
    methods
        function obj = Map2D(varargin)
        end

        function plot(obj, pos_tx, pos_rx, rot_rx_abs, area_size, aoa_act, lim_angle, show_limits)

            num_rx = size(pos_rx, 1);
            rays_abs = cell(num_rx, 1);
            % --- Plot the Tx and Rx positions
            obj.addMarkers(pos_tx(:,1), pos_tx(:,2), 'Tx', 'red');
            obj.addMarkers(pos_rx(:,1), pos_rx(:,2), 'Rx', 'blue');


            % --- For each Rx
            for rx_idx = 1:num_rx
                % Calculate the absolute rays
                rays_abs{rx_idx} = obj.calAbsRays(pos_rx(rx_idx,:), pos_tx, rot_rx_abs(rx_idx), aoa_act(rx_idx), lim_angle);
                % Plot the rays connecting the Tx and Rx
                % fplot(@(x) rays_abs{rx_idx}.doa_slope*x + rays_abs{rx_idx}.doa_shift, rays_abs{rx_idx}.lim, 'r', 'LineWidth', 2); hold on;
                % Draw connecting ray from Tx to Rx
                plot([pos_tx(1), pos_rx(rx_idx,1)], [pos_tx(2), pos_rx(rx_idx,2)], 'g-', 'LineWidth', 2);
                fplot(@(x) rays_abs{rx_idx}.centre_slope*x + rays_abs{rx_idx}.centre_shift, rays_abs{rx_idx}.lim, 'k--', 'LineWidth', 0.25); hold on;
                % add the true aoa in the plot
                text(pos_rx(rx_idx,1)+1.5, pos_rx(rx_idx,2)-5, sprintf('AoA: %.2f°', aoa_act(rx_idx)), 'Color', 'blue', 'FontSize', 12); hold on;
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
            
            % --- Plot the navigation rays
            if show_limits
                for rx_idx = 1:num_rx
                    fplot(@(x) rays_abs{rx_idx}.cw_slope*x + rays_abs{rx_idx}.cw_shift, [rays_abs{rx_idx}.lim], 'k--', 'LineWidth', 0.25); hold on;
                    fplot(@(x) rays_abs{rx_idx}.ccw_slope*x + rays_abs{rx_idx}.ccw_shift, [rays_abs{rx_idx}.lim], 'k--', 'LineWidth', 0.25); hold on;
                end
            end

            % Other format settings
            xlabel('X Position (m)'); ylabel('Y Position (m)');
            xlim([0 area_size]); ylim([0 area_size]);
            title('Map Visualisation'); grid on; hold off;
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
            if nargin == 6
                abs_ray.centre_slope = tand(rot_abs_rx);
                abs_ray.centre_shift = pos_abs_rx(2) - abs_ray.centre_slope * pos_abs_rx(1);
                abs_ray.cw_slope = tand(rot_abs_rx - rot_lim);
                abs_ray.cw_shift = pos_abs_rx(2) - abs_ray.cw_slope * pos_abs_rx(1);
                abs_ray.ccw_slope = tand(rot_abs_rx + rot_lim);
                abs_ray.ccw_shift = pos_abs_rx(2) - abs_ray.ccw_slope * pos_abs_rx(1);
                abs_ray.lim = obj.calRayPlotLim(pos_abs_rx, pos_abs_tx, abs_ray.doa_slope);
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