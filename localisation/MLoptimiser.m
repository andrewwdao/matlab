classdef MLoptimiser < handle
    % Properties to store necessary data for optimization
    properties
        l4c                 % Instance of Likelihood4Coordinates
        pos_rx             % Receiver positions
        rot_abs            % Absolute rotations
        w                  % Received signals (cell array)
        ELEMENT_NUM        % Number of ULA elements
        nPower_model       % Noise power for the model
        area_size          % Size of the search area
        grid_points        % Number of points in each direction for the coarse grid
    end
    
    % Methods of the class
    methods
        % Constructor to initialize the object
        function obj = MLoptimiser(l4c, pos_rx, rot_abs, w, ELEMENT_NUM, nPower_model, area_size, grid_points)
            obj.l4c = l4c;
            obj.pos_rx = pos_rx;
            obj.rot_abs = rot_abs;
            obj.w = w;
            obj.ELEMENT_NUM = ELEMENT_NUM;
            obj.nPower_model = nPower_model;
            obj.area_size = area_size;
            obj.grid_points = grid_points;
        end
        
        function [optCoord, L_peak] = findMaxLikelihood(obj)
            % findMaxLikelihood: Find the coordinates that maximize the likelihood
            % Outputs:
            %   optCoord - Optimal coordinates [x, y]
            %   L_peak   - Peak likelihood value
            
            % Generate a coarse grid of initial guesses
            x = linspace(0, obj.area_size, obj.grid_points);
            y = linspace(0, obj.area_size, obj.grid_points);
            [X, Y] = meshgrid(x, y);
            initial_guesses = [X(:), Y(:)]; % List of [x, y] coordinates
            
            % Initialize variables for the best result
            best_L = Inf;    % Best objective value (to minimize)
            best_coor = [0, 0]; % Best coordinate
            found_solution = false;
            
            % Optimization options (silent mode with FunValCheck)
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'FunValCheck', 'on');
            
            % Define bounds
            lb = [0, 0];
            ub = [obj.area_size, obj.area_size];
            
            % Loop over each initial guess
            for i = 1:size(initial_guesses, 1)
                coor0 = initial_guesses(i, :);
                % Define the objective function (negative likelihood for minimization)
                objective = @(coor) -abs(obj.l4c.fminconCalculateLikelihood(coor, ...
                    obj.pos_rx, obj.rot_abs, obj.w, obj.ELEMENT_NUM, obj.nPower_model));
                
                try
                    % Evaluate objective at initial point
                    test_val = objective(coor0);
                    
                    % Handle special cases
                    if isnan(test_val)
                        fprintf('Skipping initial point [%.2f, %.2f]: Objective is NaN.\n', coor0(1), coor0(2));
                        continue;
                    elseif isinf(test_val) && test_val > 0
                        fprintf('Skipping initial point [%.2f, %.2f]: Objective is +Inf.\n', coor0(1), coor0(2));
                        continue;
                    elseif isinf(test_val) && test_val < 0
                        fprintf('Found potential solution at [%.2f, %.2f]: Objective is -Inf.\n', coor0(1), coor0(2));
                        % -Inf is the best possible value for minimization
                        if -Inf < best_L
                            best_L = -Inf;
                            best_coor = coor0;
                            found_solution = true;
                        end
                        continue; % No need to run fmincon
                    end
                    
                    % If test_val is finite, run fmincon
                    [opt_coor, fval] = fmincon(objective, coor0, [], [], [], [], lb, ub, [], options);
                    fprintf('Optimization from [%.2f, %.2f] converged to [%.2f, %.2f] with value %.4f.\n', ...
                            coor0(1), coor0(2), opt_coor(1), opt_coor(2), fval);
                    
                    % Update best solution if better
                    if fval < best_L
                        best_L = fval;
                        best_coor = opt_coor;
                        found_solution = true;
                    end
                catch ME
                    fprintf('Error at initial point [%.2f, %.2f]: %s\n', coor0(1), coor0(2), ME.message);
                    continue;
                end
            end
            
            % Check if any valid solution was found
            if ~found_solution
                error('No valid initial points found. Cannot determine maximum likelihood estimate.');
            end
            
            % Return the best coordinate and the corresponding likelihood
            optCoord = best_coor;
            if isinf(best_L) && best_L < 0
                L_peak = Inf; % Objective = -L, so -Inf corresponds to L = Inf
            else
                L_peak = -best_L; % Convert back to positive likelihood
            end
        end
    end
end