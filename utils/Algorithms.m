classdef Algorithms < handle
    properties
        l4c          % Likelihood4Coordinates object
        optimiser    % Optimisers object
    end
    
    methods
        function obj = Algorithms(l4c, optimiser)
            % Constructor
            if nargin < 2
                error('Algorithms requires Likelihood4Coordinates and Optimisers objects');
            end
            obj.l4c = l4c; % Likelihood4Coordinates object
            obj.optimiser = optimiser; % Optimisers object
        end
        
        function error = errorDistance(~, pos_act, pos_est)
            % Calculate distance error between true and estimated positions
            error = norm(pos_act-pos_est); % sqrt((pos_act(1)-pos_est(1))^2 + (pos_act(2)-pos_est(2))^2);
        end

        
        function abs_ray = calRay(~, pos_rx_abs, rot_rx_abs, aoa)
            % Calculate ray parameters in absolute coordinates based on receiver position, rotation, and angle of arrival
            %
            % Inputs:
            %   pos_rx_abs - [x, y] position coordinates of the receiver in absolute reference frame
            %   rot_rx_abs - rotation angle of the receiver in absolute reference frame (degrees)
            %   aoa - angle of arrival in the receiver's local reference frame (degrees)
            %
            % Outputs:
            %   abs_ray - struct containing ray parameters in absolute coordinates
            %     .doa_slope - slope of the direction of arrival line
            %     .doa_shift - y-intercept of the direction of arrival line
            %
            % The function converts the local angle of arrival to a line equation in the absolute
            % reference frame, represented in the form: y = slope*x + shift
            abs_ray.doa_slope = tand(rot_rx_abs + aoa);
            abs_ray.doa_shift = pos_rx_abs(2) - abs_ray.doa_slope * pos_rx_abs(1);
        end

        function result = calIntersect(~, abs_ray1, abs_ray2)
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

        function [pos_est, error] = DoAintersect(obj, pos_rx, rot_abs, y_centralised, doa_estimator, pos_act)
            % DoAintersect Calculate the intersection point of DoA rays and estimate the transmitter position
            %
            % Inputs:
            %   pos_rx        - Receiver positions [2 × 2]
            %   rot_abs       - Absolute rotation angles [2 × 1]
            %   y_centralised - Cell array of received signals {2 × 1}
            %   doa_estimator - Function handle to DOA estimation algorithm
            %   pos_act       - True transmitter position [1 × 2] (optional, for error calculation)
            %
            % Outputs:
            %   pos_est       - Estimated transmitter coordinates [x, y]
            %   error         - Distance error between estimated and true position (if pos_act provided)
            
            % --- Estimate the AoA and the ray to that AoA for each receiver
            aoa_est = zeros(2, 1);           % Pre-allocate for Relative AoA estimation
            rays_abs = cell(2, 1);               % Pre-allocate for absolute rays
            for idx_rx = 1:2
                aoa_est(idx_rx, 1) = doa_estimator(y_centralised{idx_rx}).aoa_est;
                rays_abs{idx_rx, 1} = obj.calRay(pos_rx(idx_rx,:), rot_abs(idx_rx), aoa_est(idx_rx, 1));
            end

            % --- Calculate the aoa intersection point and the error distance
            pos_est = obj.calIntersect(rays_abs{1, 1}, rays_abs{2, 1});
            
            % Calculate error if true position is provided
            if nargin > 5 && ~isempty(pos_act)
                error = obj.errorDistance(pos_act, [pos_est.x, pos_est.y]);
            else
                error = [];
            end
        end

        function [pos_est, opt_val, error] = MLOptwGrid(obj, pos_rx, rot_abs, y_centralised, element_num, nPower, lb, ub, grid_density, pos_act)
            % Perform ML optimization using grid search followed by fmincon
            %
            % Inputs:
            %   pos_rx       - Receiver positions [2 × 2]
            %   rot_abs      - Absolute rotation angles [2 × 1]
            %   y_centralised - Cell array of received signals {2 × 1}
            %   element_num  - Number of array elements
            %   nPower       - Noise power
            %   pos_act     - True transmitter position [1 × 2] (optional, for error calculation)
            %   grid_bounds  - [lower_bound (list), upper_bound (list)] for x and y coordinates
            %   grid_density - Grid density for initial search
            %
            % Outputs:
            %   pos_est     - Optimized coordinates [x, y]
            %   error        - Distance error from true position (if pos_act provided)
            
            % Create objective function (convert maximization to minimization)
            obj2min = @(coor) -obj.l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, y_centralised', element_num, nPower);
            
            % Run grid-based optimization
            [pos_est, opt_val] = obj.optimiser.gridFmincon2D(obj2min, {}, ...
                lb, ub, grid_density);
            
            % Calculate error if true position is provided
            if nargin > 9 && ~isempty(pos_act)
                error = obj.errorDistance(pos_act, pos_est);
            else
                error = [];
            end
        end
        
        function [pos_est, opt_val, error] = MLOpt4mDoAtriage(obj, pos_rx, rot_abs, y_centralised, element_num, nPower, lb, ub, doa_estimator, pos_act)
            % MLOpt4mDoAtriage Maximum likelihood optimization for DOA estimate with AoA triage
            %
            % Performs maximum likelihood (ML) optimization using AoA estimates as initial point
            % to refine the transmitter position estimation. The function first estimates AoA at
            % 2 first receivers, calculates the intersection point of rays, and uses this as the 
            % initial point for ML optimization.
            %
            % Inputs:
            %   pos_rx        - Receiver positions [2 × 2], each row containing [x, y] coordinates
            %   rot_abs       - Absolute rotation angles of the receivers [2 × 1] in radians
            %   y_centralised - Cell array of received signals {2 × 1}
            %   element_num   - Number of elements in the antenna array
            %   nPower        - Noise power estimate
            %   lb            - Lower bounds for optimization [x_min, y_min]
            %   ub            - Upper bounds for optimization [x_max, y_max]
            %   doa_estimator - Function handle to DOA estimation algorithm
            %   pos_act       - True transmitter position [1 × 2] (optional, for error calculation)
            %
            % Outputs:
            %   pos_est       - Estimated transmitter coordinates [x, y]
            %   error         - Distance error between estimated and true position (if pos_act provided)
            [aoa_intersect, ~] = obj.DoAintersect(pos_rx, rot_abs, y_centralised, @(signal) doa_estimator(signal));
            % aoa_est = zeros(2, 1);           % Pre-allocate for Relative AoA estimation
            % rays_abs = cell(2, 1);               % Pre-allocate for absolute rays
            % % --- Estimate the AoA and the ray to that AoA for each receiver
            % for idx_rx = 1:2
            %     aoa_est(idx_rx, 1) = doa_estimator(y_centralised{idx_rx}).aoa_est;
            %     rays_abs{idx_rx, 1} = obj.calRay(pos_rx(idx_rx,:), rot_abs(idx_rx), aoa_est(idx_rx, 1));
            % end

            % --- Create objective function (convert maximization to minimization)
            obj2min = @(coor) -obj.l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, y_centralised', element_num, nPower);
            % --- Calculate the aoa intersection point and the error distance
            % aoa_intersect = obj.calIntersect(rays_abs{1, 1}, rays_abs{2, 1});
            [pos_est, opt_val] = obj.optimiser.findMinWithInitialPoint(obj2min, {}, lb, ub, [aoa_intersect.x, aoa_intersect.y]);
            
            % Calculate error if true position is provided
            if nargin > 9 && ~isempty(pos_act)
                error = obj.errorDistance(pos_act, pos_est);
            else
                error = [];
            end
        end

        function [centroid, pos_est, opt_val, error] = MLOpt4mCentroid(obj, pos_rx, rot_abs, y_centralised, element_num, nPower, lb, ub, doa_estimator, pos_act)
            % Maximum likelihood optimization using centroid of all ray intersections
            %
            % Uses the centroid of all DoA ray intersections as the initial point for ML optimization.
            % This approach leverages information from all available receivers, not just the first two.
            %
            % Inputs:
            %   pos_rx        - Receiver positions [N × 2], each row containing [x, y] coordinates
            %   rot_abs       - Absolute rotation angles of the receivers [N × 1] in radians
            %   y_centralised - Cell array of received signals {N × 1}
            %   element_num   - Number of elements in the antenna array
            %   nPower        - Noise power estimate
            %   lb            - Lower bounds for optimization [x_min, y_min]
            %   ub            - Upper bounds for optimization [x_max, y_max]
            %   doa_estimator - Function handle to DOA estimation algorithm
            %   pos_act       - True transmitter position [1 × 2] (optional, for error calculation)
            %
            % Outputs:
            %   pos_est       - Estimated transmitter coordinates [x, y]
            %   error         - Distance error between estimated and true position (if pos_act provided)
            
            num_rx = size(pos_rx, 1);
            aoa_est = zeros(num_rx, 1);          % Pre-allocate for Relative AoA estimation
            rays_abs = cell(num_rx, 1);          % Pre-allocate for absolute rays
            
            
            % --- Estimate the AoA and the ray to that AoA for each receiver
            for idx_rx = 1:num_rx
                aoa_est(idx_rx, 1) = doa_estimator(y_centralised{idx_rx}).aoa_est;
                rays_abs{idx_rx} = obj.calRay(pos_rx(idx_rx,:), rot_abs(idx_rx), aoa_est(idx_rx, 1));
            end
            
            % --- Calculate all intersection points between rays
            intersection_points = struct('x', {}, 'y', {});
            valid_count = 0;
            
            for i = 1:(num_rx-1)
                for j = (i+1):num_rx
                    % Skip if ray slopes are nearly identical (parallel lines)
                    if abs(rays_abs{i}.doa_slope - rays_abs{j}.doa_slope) < 1e-10
                        continue;
                    end
                    
                    % Calculate intersection
                    intersection = obj.calIntersect(rays_abs{i}, rays_abs{j});
                    
                    % Check if intersection point is within bounds
                    if intersection.x >= lb(1) && intersection.x <= ub(1) && ...
                       intersection.y >= lb(2) && intersection.y <= ub(2)
                        valid_count = valid_count + 1;
                        intersection_points(valid_count) = intersection;
                    end
                end
            end
            
            % --- Calculate centroid of all valid intersection points
            if valid_count > 0
                x_sum = 0;
                y_sum = 0;
                for i = 1:valid_count
                    x_sum = x_sum + intersection_points(i).x;
                    y_sum = y_sum + intersection_points(i).y;
                end
                centroid = [x_sum/valid_count, y_sum/valid_count];
            else
                % Fallback if no valid intersections found: use midpoint of area
                centroid = [(lb(1) + ub(1))/2, (lb(2) + ub(2))/2];
            end

            % Create objective function (convert maximization to minimization)
            obj2min = @(coor) -obj.l4c.likelihoodFromCoorSet(coor, pos_rx, rot_abs, y_centralised', element_num, nPower);
            % --- Use centroid as initial point for ML optimization
            [pos_est, opt_val] = obj.optimiser.findMinWithInitialPoint(obj2min, {}, lb, ub, centroid);
            
            % Calculate error if true position is provided
            if nargin > 9 && ~isempty(pos_act)
                error = obj.errorDistance(pos_act, pos_est);
            else
                error = [];
            end
        end
    end
end