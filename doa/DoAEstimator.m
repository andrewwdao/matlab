classdef DoAEstimator < handle
    properties
        received_signal % Received signal at the sensor array
        tx_num % Number of transmitters
        lambda % Wavelength (m)
        element_num % Number of elements in the ULA
        element_spacing % Distance between antenna elements (normalised to lambda units)
        sweeping_angle % Angle range for sweeping to find the AoA
        steer_vect % Steering vector that simulate the effects of the array on the received signal
        aoa_act % Actual Angle of Arrival
    end

    methods
        function obj = DoAEstimator(received_signal, tx_num, lambda, element_num, element_spacing, sweeping_angle, aoa_act)
            % --- Constructor to initialize the properties
            % steer_vect: complex Nx1 (Number of elements)x1
            obj.received_signal = received_signal;  % Received signal at the sensor array
            obj.tx_num = tx_num;  % Number of transmitters
            obj.lambda = lambda;  % Wavelength (m)
            obj.element_num = element_num;  % Number of elements in the ULA
            obj.element_spacing = element_spacing;  % Distance between antenna elements (normalised to lambda units))
            obj.sweeping_angle = sweeping_angle;  % Angle range for sweeping to find the AoA
            obj.aoa_act = aoa_act;
            obj.steer_vect = zeros(element_num, length(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                obj.steer_vect(:, i) = exp(-2j * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.sweeping_angle(i)) / obj.lambda);
            end
        end

        function logger(obj, purpose, values)
            disp('--------------------------------------------------------');
            disp('Estimated Angles of Arrival');
            disp(purpose);
            disp('--------------------------------------------------------');
            disp(array2table(...
                values, ...% table data
                'RowNames', cellstr(strcat('RX', num2str((1:1)'))), ...
                'VariableNames', cellstr(strcat('TX', num2str((1:obj.tx_num)')))));
        end

        function square_err = calculate_square_error(obj, aoa_est)
            % Calculate the square error between the estimated and true AoA
            square_err = (aoa_est - obj.aoa_act).^2;
        end

        function result = parse_output(obj, type, spectrum)
            % Convert the spectrum to dB scale
            result.spectrum_dB = 10 * log10(abs(spectrum));
            % Find the peaks in the spectrum
            % [~, peak_indices] = findpeaks(spectrum_dB,'SortStr','descend');
            % aoa_est = obj.sweeping_angle(peak_indices(1:obj.tx_num));
            [~, max_idx] = max(result.spectrum_dB);
            result.aoa_est = obj.sweeping_angle(max_idx);
            result.square_err = obj.calculate_square_error(result.aoa_est);
            % obj.logger(type, aoa_est);
        end


        function result = ML_sync(obj, s_t)
            % Assume received signal: synchoronous signal with noise
            % Model: y = steering_vector*s_t + noise;
            % (know the transmitted symbols before channel and array effects)
            % s_t: complex 1xT (1 TX)x(Number of samples)
            % obj.received_signal: complex NxT (Number of elements)x(Number of samples)
            % steering_vector: complex Nx1 (Number of elements)x1
            spectrum = zeros(size(obj.sweeping_angle));
            steer_vect_local = obj.steer_vect;
            received_signal_local = obj.received_signal;
            parfor i = 1:length(obj.sweeping_angle)
                for t = 1:size(s_t,2) % consider each time instance separately
                    spectrum(i) = spectrum(i)+ real(received_signal_local(:,t)' * steer_vect_local(:, i) * s_t(:,t)); %#ok<PFBNS>
                end
            end
            % Parse the output to readable format
            result = obj.parse_output('Sync ML', spectrum);
        end
        
        function result = ML_async(obj, s_t)
            % Assume received signal: asynchoronous signal with noise
            % Model: y = steering_vector*exp(j\theta) + noise;
            % (know the transmitted symbols before channel and array effects)
            % s_t: complex 1xT (1 TX)x(Number of samples)
            % obj.received_signal: complex NxT (Number of elements)x(Number of samples)
            % steering_vector: complex Nx1 (Number of elements)x1
            spectrum = zeros(size(obj.sweeping_angle));
            steer_vect_local = obj.steer_vect;
            received_signal_local = obj.received_signal;
            parfor i = 1:length(obj.sweeping_angle)
                for t = 1:size(s_t,2) % consider each time instance separately
                    spectrum(i) = spectrum(i) + received_signal_local(:,t)' * steer_vect_local(:, i) * steer_vect_local(:, i)' * received_signal_local(:,t); %#ok<PFBNS>
                end
            end
            % Parse the output to readable format
            result = obj.parse_output('ASync ML', spectrum);
        end

        function result = BF(obj)
            % --- Conventional Beamforming Algorithm
            % Similar to the ML estimator for asynchoronous received signal:
            % y = steering_vector*exp(j\theta) + noise;
            % --- Calculate the covariance matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2); % y*y^H/N
            spectrum = zeros(size(obj.sweeping_angle));
            steer_vect_local = obj.steer_vect;
            parfor i = 1:length(obj.sweeping_angle)
                % Form MVDR denominator matrix from noise subspace eigenvectors
                spectrum(i) = steer_vect_local(:, i)' * R * steer_vect_local(:, i);  %steering_vector' * inv(R) * steering_vector;
            end
            % Parse the output to readable format
            result = obj.parse_output('BF', spectrum);
        end
        
        function result = MVDR(obj)
            % --- MVDR Algorithm - or Minimum Variance Distortionless Response - Maximum Likelihood
            % --- Calculate the covariance matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2); % y*y^H/N
            % Add regularization to the covariance matrix (diagonal loading)
            if rcond(R) < 1e-15  % if R is a (near)singular matrix (non invertible)
                R = R + eye(size(R));
            end
            spectrum = zeros(size(obj.sweeping_angle));
            steer_vect_local = obj.steer_vect;
            parfor i = 1:length(obj.sweeping_angle)
                % Form MVDR denominator matrix from noise subspace eigenvectors
                denom = (steer_vect_local(:, i)' / R) * steer_vect_local(:, i); % steering_vector' * inv(R) * steering_vector;
                spectrum(i) = 1 ./ denom;
            end
            % Parse the output to readable format
            result = obj.parse_output('MVDR', spectrum);
        end

        function result = MUSIC(obj)
            % --- MUSIC Algorithm
            % --- Calculate the covariance matrix - sample correlation matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2); % y*y^H/N
            % Perform eigenvalue decomposition
            [eigenvectors, eigenvalues] = eig(R);
            % Sort eigenvalues and eigenvectors
            [~, idx] = sort(diag(eigenvalues), 'descend');
            eigenvectors = eigenvectors(:, idx); % short the largest eigenvectors to the largest eigenvalues first
            % Determine the noise subspace
            noise_subspace = eigenvectors(:, obj.tx_num+1:end);
            % Compute the MUSIC spectrum
            spectrum = zeros(size(obj.sweeping_angle));
            steer_vect_local = obj.steer_vect;
            parfor i = 1:length(obj.sweeping_angle)
                % Form MUSIC denominator matrix from noise subspace eigenvectors
                % another implementation for the denominator:
                %(steering_vector' * (noise_subspace * noise_subspace') * steering_vector +eps(1));
                denom = sum(abs(noise_subspace' * steer_vect_local(:, i)).^2, 1);
                denom = denom+eps(1); % add a small positive constant to prevent division by zero. 9.44 in [1]
                spectrum(i) = 1 ./ denom;
            end
            % Parse the output to readable format
            result = obj.parse_output('MUSIC', spectrum);
        end

        function result = MUSIC_op(obj)
            % --- MUSIC Algorithm
            % Calculate covariance matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2);
            
            % Eigenvalue decomposition
            [eigenvectors, eigenvalues] = eig(R);
            [~, idx] = sort(diag(eigenvalues), 'descend');
            eigenvectors = eigenvectors(:, idx);
            noise_subspace = eigenvectors(:, obj.tx_num+1:end);
            
            % Array parameters
            Nr = size(obj.received_signal, 1); % Number of elements
            d = obj.element_spacing; % Element spacing (normalized by wavelength)
            
            % Step 1: Coarse grid search for initial guesses
            coarse_angles = linspace(-90, 90, 37); % 5-degree resolution
            coarse_spectrum = zeros(size(coarse_angles));
            for i = 1:length(coarse_angles)
                theta = coarse_angles(i);
                steer_vec = exp(-1i * 2 * pi * d * (0:Nr-1)' * sind(theta));
                denom = sum(abs(noise_subspace' * steer_vec).^2) + eps(1);
                coarse_spectrum(i) = 1 / denom;
            end
            
            % Find peaks in coarse spectrum
            [~, peak_locs] = findpeaks(coarse_spectrum, 'SortStr', 'descend', 'NPeaks', obj.tx_num);
            initial_guesses = coarse_angles(peak_locs);
            
            % Step 2: Refine angles using optimization
            optimized_angles = zeros(1, obj.tx_num);
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');
            
            % Objective function (minimize negative MUSIC spectrum)
            music_obj = @(theta) -1/(sum(abs(noise_subspace' * ...
                                  exp(-1i * 2 * pi * d * (0:Nr-1)' * sind(theta))).^2) + eps(1));
            
            % Refine each initial guess
            for k = 1:obj.tx_num
                [theta_opt, ~] = fmincon(music_obj, initial_guesses(k), [], [], [], [], -90, 90, [], options);
                optimized_angles(k) = theta_opt;
            end
            
            % Optional: Compute full spectrum for visualization
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                denom = sum(abs(noise_subspace' * obj.steer_vect(:, i)).^2) + eps(1);
                spectrum(i) = 1 / denom;
            end
            
            % Parse output with optimized angles
            result = obj.parse_output('MUSIC', spectrum);
        end
    end
end
