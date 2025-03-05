classdef gridOptimiser < handle
    % Properties to store necessary data for optimization
    properties
    end
    
    methods
        % Constructor to initialize the object
        function obj = gridOptimiser(varargin)
        end
        
        % Method to find the min of the objective function using fmincon
        function [opt_var, opt_val] = fmincon1D(~, objective_to_minimise, extra_args, lb, ub, grid_points)
            % Generate a coarse grid of initial guesses within bounds
            initial_guesses = linspace(lb, ub, grid_points);
            
            % Initialize variables for the best result
            opt_val = Inf;  % Best objective value (to minimize)
            opt_var = [];   % Best coordinate
            found_solution = false;
            
            % Optimization options (silent mode with FunValCheck)
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'FunValCheck', 'on');
            
            % Loop over each initial guess
            for i = 1:length(initial_guesses)
                angle0 = initial_guesses(i);
                % Define the objective function for fmincon
                objective_for_fmincon = @(angle) objective_to_minimise(angle, extra_args{:});
                
                try
                    % Evaluate objective at initial point
                    test_val = objective_for_fmincon(angle0);
                    
                    % Handle special cases
                    if isnan(test_val) || isinf(test_val) && test_val > 0
                        fprintf('Skipping initial point [%.2f, %.2f]: Objective is NaN/+Inf.\n', angle0);
                        continue;
                    elseif isinf(test_val) && test_val < 0
                        fprintf('Found potential solution at [%.2f, %.2f]: Objective is -Inf.\n', angle0);
                        % -Inf is the best possible value for minimization
                        opt_val = -Inf;
                        opt_var = angle0;
                        found_solution = true;
                        continue; % No need to run fmincon
                    end
                    
                    % If test_val is finite, run fmincon
                    [opt_coor, fval] = fmincon(objective_for_fmincon, angle0, [], [], [], [], lb, ub, [], options);
                    % fprintf('Initial point %.2f converged to %.2f with value %.4f.\n', angle0, opt_coor, fval);
                    
                    % Update best solution if better
                    if fval < opt_val
                        opt_val = fval;
                        opt_var = opt_coor;
                        found_solution = true;
                    end
                catch ME
                    fprintf('Error at initial point %.2f: %s\n', angle0, ME.message);
                    continue;
                end
            end
            
            % Check if any valid solution was found
            if ~found_solution
                error('No valid initial points found. Cannot determine local minimum value.');
            end
        end

        % Method to find the min of the objective function of 2 values using fmincon
        function [optCoord, L_peak] = fmincon2D(~, objective_to_minimise, extra_args, lb, ub, grid_points)
            % Generate a coarse grid of initial guesses within bounds
            x = linspace(lb(1), ub(1), grid_points);
            y = linspace(lb(2), ub(2), grid_points);
            [X, Y] = meshgrid(x, y);
            initial_guesses = [X(:), Y(:)]; % List of [x, y] coordinates
            
            % Initialize variables for the best result
            best_fval = Inf;  % Best objective value (to minimize)
            best_coor = [];   % Best coordinate
            found_solution = false;
            
            % Optimization options (silent mode with FunValCheck)
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'FunValCheck', 'on');
            
            % Loop over each initial guess
            for i = 1:size(initial_guesses, 1)
                coor0 = initial_guesses(i, :);
                % Define the objective function for fmincon
                objective_for_fmincon = @(coor) objective_to_minimise(coor, extra_args{:});
                
                try
                    % Evaluate objective at initial point
                    test_val = objective_for_fmincon(coor0);
                    
                    % Handle special cases
                    if isnan(test_val) || isinf(test_val) && test_val > 0
                        fprintf('Skipping initial point [%.2f, %.2f]: Objective is NaN/+Inf.\n', coor0(1), coor0(2));
                        continue;
                    elseif isinf(test_val) && test_val < 0
                        fprintf('Found potential solution at [%.2f, %.2f]: Objective is -Inf.\n', coor0(1), coor0(2));
                        % -Inf is the best possible value for minimization
                        best_fval = -Inf;
                        best_coor = coor0;
                        found_solution = true;
                        continue; % No need to run fmincon
                    end
                    
                    % If test_val is finite, run fmincon
                    [opt_coor, fval] = fmincon(objective_for_fmincon, coor0, [], [], [], [], lb, ub, [], options);
                    % fprintf('Initial point [%.2f, %.2f] converged to [%.2f, %.2f] with value %.4f.\n', ...
                            % coor0(1), coor0(2), opt_coor(1), opt_coor(2), fval);
                    
                    % Update best solution if better
                    if fval < best_fval
                        best_fval = fval;
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
                error('No valid initial points found. Cannot determine local minimum value.');
            end
            
            % Return the best coordinate and the corresponding maximized value
            optCoord = best_coor;
            L_peak = -best_fval;  % Since objective_for_fmincon = -objective_to_minimise
        end
    end
end
