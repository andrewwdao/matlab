
classdef Map2D < handle
    properties
        POS_RX = [
            0, 0; % 1 - Bottom-left corner
            100, 0; % 2 - Bottom-right corner
            100, 100; % 3 - Top-right corner
            0,100; % 4 - Top-left corner
            0, 50; % 5 - Left edge
            50, 0; % 6 - Bottom edge
            100, 50; % 7 - Right edge
            50, 100; % 8 - Top edge
            25,100; % 9 - Top edge
            0, 75; % 10 - Left edge
            0, 25; % 11 - Left edge
            25,0; % 12 - Bottom edge
            75,0; % 13 - Bottom edge
            100, 25; % 14 - Right edge
            100, 75; % 15 - Right edge
        ];
        POS_TX = [
            50, 50; % 1 - Center of the area
            25, 25; % 2 - Bottom-left corner
            75, 25; % 3 - Bottom-right corner
            75, 75; % 4 - Top-right corner
            25, 75; % 5 - Top-left corner
        ];

    end

    methods
        function obj = Map2D(varargin)
            if nargin == 3
                obj.POS_RX = obj.genSquarePerimeterNodes(varargin{1}, varargin{2}, varargin{3});
            end
        end

        function POS_RX = genSquarePerimeterNodes(~, start_point, end_point, num_positions)
            % Generate positions in a zigzag pattern around a square perimeter
            % Points are placed at corners first, then midpoints of edges, then recursively subdivide
            
            % Extract square coordinates
            x1 = start_point(1);
            y1 = start_point(2);
            x2 = end_point(1);
            y2 = end_point(2);
            
            % Initialize with empty array
            POS_RX = zeros(num_positions, 2);
            
            % Start with the four corners in specific order
            corners = [
                x1, y1;     % Bottom-left (start point)
                x2, y1;     % Bottom-right
                x2, y2;     % Top-right
                x1, y2      % Top-left
            ];
            
            % Add the corners as the first points
            num_corners = size(corners, 1);
            POS_RX(1:num_corners, :) = corners;
            
            % Define the four edges (in order: left, bottom, right, top)
            edges = [
                [x1, y1; x1, y2];  % Left edge
                [x1, y1; x2, y1];  % Bottom edge
                [x2, y1; x2, y2];  % Right edge
                [x1, y2; x2, y2]   % Top edge
            ];
            
            % Estimate maximum number of edges needed (conservative estimate)
            max_edges = 2 * (num_positions - num_corners) + 4;
            
            % Queue to store edges to be subdivided
            edge_queue = cell(1, max_edges);
            queue_size = 4;
            for i = 1:4
                edge_queue{i} = edges(i*2-1:i*2, :);
            end
            
            % Current point index
            point_idx = num_corners + 1;
            
            % Process each edge in the queue until we have enough points
            current_edge = 1;
            while point_idx <= num_positions && current_edge <= queue_size
                % Get current edge
                edge = edge_queue{current_edge};
                
                % Calculate midpoint
                midpoint = [(edge(1,1) + edge(2,1))/2, (edge(1,2) + edge(2,2))/2];
                
                % Add midpoint to the positions
                POS_RX(point_idx, :) = midpoint;
                point_idx = point_idx + 1;
                
                % Split current edge into two new edges and add to queue
                if point_idx <= num_positions
                    queue_size = queue_size + 1;
                    edge_queue{queue_size} = [edge(1,:); midpoint];
                    queue_size = queue_size + 1;
                    edge_queue{queue_size} = [midpoint; edge(2,:)];
                end
                
                % Move to next edge in queue
                current_edge = current_edge + 1;
            end
            
            % Trim to requested number of positions
            POS_RX = POS_RX(1:num_positions, :);
        end

        function plots(obj, pos_tx, pos_rx, rot_rx_abs, area_size, aoa_act, lim_angle, flags, aoa_est_cell)
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
                    plot([pos_tx(tx_idx,1), pos_rx(rx_idx,1)], [pos_tx(tx_idx,2), pos_rx(rx_idx,2)], 'Color', [0.95, 0.65, 0.3], 'LineWidth', 1); hold on;
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
                    text(pos_rx(rx_idx,1)+1.5, pos_rx(rx_idx,2)-4.5, sprintf('AoA: %.2f°', aoa_act(rx_idx)), 'Color', 'blue', 'FontSize', 12); hold on;
                    % add the estimated aoa in the plot if available
                    if nargin == 9
                        text(pos_rx(rx_idx,1)+1.5, pos_rx(rx_idx,2)-12.5, sprintf('Est. AoA: %.2f°', aoa_est_cell{1, rx_idx}), 'Color', 'blue', 'FontSize', 12); hold on;
                    end              
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
            obj.plots(pos_tx, pos_rx, rot_rx_abs, area_size, aoa_act, lim_angle, flags, aoa_est_cell);hold on;
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

        function abs_angle = calAbsAngle(~, pos_tx, pos_rx, aoa_act)
            % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants.
            NUM_RX = size(pos_rx, 1);
            angle_rx_tx_abs = zeros(NUM_RX, 1);
            for i = 1:NUM_RX
                angle_rx_tx_abs(i) = atan2d(pos_tx(2)-pos_rx(i,2), pos_tx(1)-pos_rx(i,1));
            end
            abs_angle = angle_rx_tx_abs - aoa_act; % Absolute rotation of the receiver in degrees
        end
        
        function abs_ray = calAbsRays(obj, pos_rx_abs, pos_tx_abs, rot_rx_abs, aoa_act, rot_lim)
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
            %    pos_rx_abs - A 1x2 vector representing the absolute position of the Rx in the 2D plane [x, y].
            %    pos_tx_abs - A 1x2 vector representing the absolute position of the Tx in the 2D plane [x, y].
            %    rot_rx_abs - A scalar representing the absolute rotation of the Rx in degrees.
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
            %    pos_rx_abs = [100, 200];
            %    pos_tx_abs = [300, 400];
            %    rot_rx_abs = 45;
            %    aoa_act = 30;
            abs_ray.doa_slope = tand(rot_rx_abs + aoa_act);
            abs_ray.doa_shift = pos_rx_abs(2) - abs_ray.doa_slope * pos_rx_abs(1);
            abs_ray.lim = obj.calRayPlotLim(pos_rx_abs, pos_tx_abs, abs_ray.doa_slope);
            if nargin == 6
                abs_ray.centre_slope = tand(rot_rx_abs);
                abs_ray.centre_shift = pos_rx_abs(2) - abs_ray.centre_slope * pos_rx_abs(1);
                abs_ray.cw_slope = tand(rot_rx_abs - rot_lim);
                abs_ray.cw_shift = pos_rx_abs(2) - abs_ray.cw_slope * pos_rx_abs(1);
                abs_ray.ccw_slope = tand(rot_rx_abs + rot_lim);
                abs_ray.ccw_shift = pos_rx_abs(2) - abs_ray.ccw_slope * pos_rx_abs(1);
            end
        end

        function ray_plot_lim = calRayPlotLim(~, pos_rx_abs, pos_tx_abs, slope)
            % Calculate the limits of the ray plot based on the quadrant in which the Rx and Tx are located.
            % 
            % Inputs:
            %    pos_rx_abs - A 1x2 vector representing the absolute position of the Rx in the 2D plane [x, y].
            %    pos_tx_abs - A 1x2 vector representing the absolute position of the Tx in the 2D plane [x, y].
            %    slope      - A scalar representing the slope of the line connecting the Tx and Rx.
            % 
            % Outputs:
            %    ray_plot_lim - A 1x2 vector representing the limits of the ray plot in the 2D plane.
            % 
            % Example:
            %    pos_rx_abs = [100, 200];
            %    pos_tx_abs = [300, 400];
            %    slope = 0.5;
            if (pos_rx_abs(1) < pos_tx_abs(1) && pos_rx_abs(2) < pos_tx_abs(2)) || (pos_rx_abs(1) > pos_tx_abs(1) && pos_rx_abs(2) < pos_tx_abs(2))
                lim_range = [0, 10000]; % Third quadrant or first quadrant
            else
                lim_range = [-10000, 0]; % Second quadrant or fourth quadrant
            end
            ray_plot_lim = sort(slope * lim_range + pos_rx_abs(1), "ascend"); % Absolute ray limit from a relative x limit from 0 to 100
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

        function [pos_rx, aoa_act, rot_abs] = genRXPos(obj, area_size, pos_tx, rx_num, rx_randomised, safety_dist, angle_limit, resolution)
            % Generate random positions for the Rx within the area size and a safe distance from the Tx.
            % 
            % Inputs:
            %    area_size - The size of the area in which the Rx can be placed.
            %    pos_tx - The position of the Tx in the 2D plane.
            %    safety_dist - The minimum distance between the Tx and Rx.
            %    rx_num - The number of Rx to generate.
            %    angle_limit - The current angle limit for the AoA.
            %    resolution - The resolution of the AoA.
            %
            % Outputs:
            %    pos_rx - The generated positions of the Rx in the 2D plane.
            %    aoa_act - The generated true AoA within the current angle limit.
            %    rot_abs - The absolute rotation of the receiver in degrees.
            %
            % Example:
            %    area_size = 100;
            %    pos_tx = [50, 50];
            %    safety_dist = 10;
            %    rx_num = 2;
            %    angle_limit = 30 degree;
            %    resolution = 1 degree;
            pos_rx = zeros(rx_num, 2);
            aoa_act = zeros(rx_num, 1);
            if rx_randomised
                for i = 1:rx_num
                    valid_position = false;
                    while ~valid_position
                        % Generate random position
                        pos_rx(i,:) = area_size * rand(1, 2);
                        % Check if it's far enough from TX
                        if sqrt(sum((pos_tx - pos_rx(i,:)).^2)) >= safety_dist
                            valid_position = true;
                        end
                    end
                    % Generate random true Angle of Arrival within the current angle limit
                    aoa_act(i) = -angle_limit + resolution * randi([0, 2*angle_limit/resolution], 1, 1);
                  
                end
            else % Fixed RX location with 0 deg AoA based on the number of rx
                for i = 1:rx_num
                    pos_rx(i,:) = obj.POS_RX(i,:);
                    aoa_act(i) = 0;
                end
            end
            % Calculate the absolute angle of the receiver to the transmitter with 4 quadrants
            rot_abs = obj.calAbsAngle(pos_tx, pos_rx, aoa_act);
        end

        function pos_tx = genTXPos(~, area_size, tx_num, tx_randomised)
            % Generate random positions for the Tx within the area size and a safe distance from the Rx.
            % 
            % Inputs:
            %    area_size - The size of the area in which the Tx can be placed.
            %    pos_tx - The position of the Tx in the 2D plane.
            %    safety_dist - The minimum distance between the Tx and Rx.
            %    tx_num - The number of Tx to generate.
            %
            % Outputs:
            %    pos_tx - The generated positions of the Tx in the 2D plane.
            %
            % Example:
            %    area_size = 100;
            %    pos_tx = [50, 50];
            %    safety_dist = 10;
            %    tx_num = 2;
            pos_tx = zeros(tx_num, 2);
            if tx_randomised
                for i = 1:tx_num
                    % valid_position = false;
                    % while ~valid_position
                        % Generate random position
                        pos_tx(i,:) = area_size * rand(1, 2);
                        % Check if it's far enough from RX
                        % if sqrt(sum((pos_tx - pos_tx(i,:)).^2)) >= safety_dist
                        %     valid_position = true;
                        % end
                    % end
                end
            else
                for i = 1:tx_num
                    pos_tx(i,:) = obj.POS_TX(i,:);
                end
            end
        end
    end
end