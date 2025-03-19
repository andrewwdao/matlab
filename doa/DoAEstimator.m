classdef DoAEstimator < handle
    properties
        antenna_array       % Instance of AntennaArray (e.g., ULA)
        sweeping_angle      % Angle range for sweeping to find the AoA
        steer_vect          % Precomputed steering vectors for sweeping angles
        estimation_mode = 'sweep'  % Default mode: 'sweep', 'opt', or 'both'
        grid_points = 10   % Number of grid points for optimization
        aoa_act=0             % Actual AoA for error calculation
    end
    
    methods
        % Constructor: Initialize with antenna array, sweeping angles, and optional estimation_mode
        function obj = DoAEstimator(antenna_array, sweeping_angle, aoa_act, estimation_mode, grid_points)
            % Inputs:
            %   antenna_array: Instance of AntennaArray (e.g., ULA)
            %   sweeping_angle: Array of angles for coarse grid search
            %   estimation_mode (optional): 'sweep', 'opt', or 'both'
            obj.antenna_array = antenna_array;
            obj.sweeping_angle = sweeping_angle;
            if nargin >= 3
                obj.aoa_act = aoa_act;
            end
            if nargin >= 4
                obj.estimation_mode = estimation_mode;
            end
            if nargin >= 5
                obj.grid_points = grid_points;
            end
            if ~strcmp(obj.estimation_mode, 'opt')
                % Precompute steering vectors for sweeping angles
                obj.steer_vect = antenna_array.getSteeringVector(sweeping_angle(:));
            end
        end
        
        % Utility method to log estimated angles
        function logger(~, purpose, values)
            tx_num = length(values);
            disp('--------------------------------------------------------');
            disp('Estimated Angles of Arrival');
            disp(purpose);
            disp('--------------------------------------------------------');
            disp(array2table(values, ...
                'RowNames', cellstr(strcat('RX', num2str((1:1)'))), ...
                'VariableNames', cellstr(strcat('TX', num2str((1:tx_num)')))));
        end
        
        function result = applySweeping(obj, objective)
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                spectrum(i) = objective(obj.sweeping_angle(i));
            end
            % Find peak in spectrum
            [~, max_idx] = max(spectrum);
            % Store results
            result.aoa_est = obj.sweeping_angle(max_idx);
            result.square_err = (result.aoa_est - obj.aoa_act).^2;
            result.spectrum_dB = 10 * log10(abs(spectrum)); % Convert to dB for plotting
        end

        function result = applyOptimization(obj, objective)
            % Instantiate Optimisers
            optimiser = Optimisers();
            % Objective function needs to be reversed for minimization
            objective = @(angle) -objective(angle);
            % Global search over full angle range
            lb = -90; % Lower bound angle
            ub = 90; % Upper bound angle
            % Perform optimization using Optimisers.fmincon
            [opt_angle, ~] = optimiser.gridFmincon1D(objective, {}, lb, ub, obj.grid_points);
            % Store results
            result.aoa_est = opt_angle;
            result.square_err = (opt_angle - obj.aoa_act).^2;
        end

        % Parse output and compute AoA based on estimation_mode
        function result = parse_output(obj, objective)
            if strcmp(obj.estimation_mode, 'sweep')
                result= obj.applySweeping(objective);
            elseif strcmp(obj.estimation_mode, 'opt')
                result = obj.applyOptimization(objective);
            end
        end


        function result = ML_sync(obj, received_signal, transmitted_signal)
            % Assume received signal: synchoronous signal with noise
            % Model: y = steering_vector*transmitted_signal + noise;
            % (know the transmitted symbols before channel and array effects)
            %   transmitted_signal: complex 1xT (1 TX)x(Number of samples)
            %   received_signal: complex NxT (Number of elements)x(Number of samples)
            %   steer_vec: complex Nx1 (Number of elements)x1
            steer_vec = @(theta) obj.antenna_array.getSteeringVector(theta);
            objective_to_maximise = @(theta) real(sum((received_signal' * steer_vec(theta)) .* transmitted_signal.'));
            result = obj.parse_output(objective_to_maximise);
        end
        
        function result = ML_async(obj, received_signal)
            % Assume received signal: asynchoronous signal with noise
            % Model: y = steering_vector*exp(j\theta) + noise;
            %   received_signal: complex NxT (Number of elements)x(Number of samples)
            %   steer_vec: complex Nx1 (Number of elements)x1
            steer_vec = @(theta) obj.antenna_array.getSteeringVector(theta);
            objective_to_maximise = @(theta) real(sum(received_signal' * steer_vec(theta) * steer_vec(theta)' * received_signal));
            result = obj.parse_output(objective_to_maximise);
        end

        function result = BF(obj, received_signal)
            % --- Conventional Beamforming Algorithm
            % Similar to the ML estimator for asynchoronous received signal:
            % y = steering_vector*exp(j\theta) + noise;
            steer_vec = @(theta) obj.antenna_array.getSteeringVector(theta);
            % Sample covariance matrix
            R = received_signal * received_signal' / size(received_signal, 2); % y*y^H/N
            objective_to_maximise = @(theta) steer_vec(theta)' * R * steer_vec(theta);
            result = obj.parse_output(objective_to_maximise);
        end
        
        function result = MVDR(obj, received_signal)
            % --- MVDR Algorithm - or Minimum Variance Distortionless Response - Maximum Likelihood
            steer_vec = @(theta) obj.antenna_array.getSteeringVector(theta);
            % Sample covariance matrix
            R = received_signal * received_signal' / size(received_signal, 2); % y*y^H/N
            % Add regularization to the covariance matrix (diagonal loading)
            if rcond(R) < 1e-15  % if R is a (near)singular matrix (non invertible)
                R = R + eye(size(R));
            end
            denom = @(theta) steer_vec(theta)' / R * steer_vec(theta);
            objective_to_maximise = @(theta) 1 ./ denom(theta);
            result = obj.parse_output(objective_to_maximise);
        end

        function result = MUSIC(obj, received_signal, tx_num)
            % --- MUSIC Algorithm
            steer_vec = @(theta) obj.antenna_array.getSteeringVector(theta);
            % Sample covariance matrix
            R = received_signal * received_signal' / size(received_signal, 2); % y*y^H/N
            % Perform eigenvalue decomposition
            [eigenvectors, eigenvalues] = eig(R);
            % Sort eigenvalues and eigenvectors
            [~, idx] = sort(diag(eigenvalues), 'descend');
            eigenvectors = eigenvectors(:, idx); % short the largest eigenvectors to the largest eigenvalues first
            % Determine the noise subspace
            noise_subspace = eigenvectors(:, tx_num+1:end);
            % Compute the MUSIC spectrum
            denom = @(theta) sum(abs(noise_subspace' * steer_vec(theta)).^2, 1)+eps(1); % add a small positive constant to prevent division by zero. 9.44 in [1]
            objective_to_maximise = @(theta) 1 ./ denom(theta);
            result = obj.parse_output(objective_to_maximise);
        end

        % function result = MUSIC_op(obj)
        %     % --- MUSIC Algorithm
        %     % Calculate covariance matrix
        %     R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2);
            
        %     % Eigenvalue decomposition
        %     [eigenvectors, eigenvalues] = eig(R);
        %     [~, idx] = sort(diag(eigenvalues), 'descend');
        %     eigenvectors = eigenvectors(:, idx);
        %     noise_subspace = eigenvectors(:, obj.tx_num+1:end);
            
        %     % Array parameters
        %     Nr = size(obj.received_signal, 1); % Number of elements
        %     d = obj.element_spacing; % Element spacing (normalized by wavelength)
            
        %     % Step 1: Coarse grid search for initial guesses
        %     coarse_angles = linspace(-90, 90, 37); % 5-degree resolution
        %     coarse_spectrum = zeros(size(coarse_angles));
        %     for i = 1:length(coarse_angles)
        %         theta = coarse_angles(i);
        %         steer_vec = exp(-1i * 2 * pi * d * (0:Nr-1)' * sind(theta));
        %         denom = sum(abs(noise_subspace' * steer_vec).^2) + eps(1);
        %         coarse_spectrum(i) = 1 / denom;
        %     end
            
        %     % Find peaks in coarse spectrum
        %     [~, peak_locs] = findpeaks(coarse_spectrum, 'SortStr', 'descend', 'NPeaks', obj.tx_num);
        %     initial_guesses = coarse_angles(peak_locs);
            
        %     % Step 2: Refine angles using optimization
        %     optimized_angles = zeros(1, obj.tx_num);
        %     options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');
            
        %     % Objective function (minimize negative MUSIC spectrum)
        %     music_obj = @(theta) -1/(sum(abs(noise_subspace' * ...
        %                           exp(-1i * 2 * pi * d * (0:Nr-1)' * sind(theta))).^2) + eps(1));
            
        %     % Refine each initial guess
        %     for k = 1:obj.tx_num
        %         [theta_opt, ~] = fmincon(music_obj, initial_guesses(k), [], [], [], [], -90, 90, [], options);
        %         optimized_angles(k) = theta_opt;
        %     end
            
        %     % Optional: Compute full spectrum for visualization
        %     spectrum = zeros(size(obj.sweeping_angle));
        %     for i = 1:length(obj.sweeping_angle)
        %         denom = sum(abs(noise_subspace' * obj.steer_vect(:, i)).^2) + eps(1);
        %         spectrum(i) = 1 / denom;
        %     end
            
        %     % Parse output with optimized angles
        %     result = obj.parse_output('MUSIC', spectrum);
        % end
    end
end
