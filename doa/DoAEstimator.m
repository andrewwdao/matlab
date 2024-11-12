classdef DoAEstimator < handle
    properties
        received_signal % Received signal at the sensor array
        tx_num % Number of transmitters
        lambda % Wavelength (m)
        element_num % Number of elements in the ULA
        element_spacing % Distance between antenna elements (normalised to lambda units)
        sweeping_angle % Angle range for sweeping to find the AoA
        steer_vect % Steering vector that simulate the effects of the array on the received signal
        act_aoa % Actual Angle of Arrival
    end

    methods
        function obj = DoAEstimator(received_signal, tx_num, lambda, element_num, element_spacing, sweeping_angle, act_aoa)
            % --- Constructor to initialize the properties
            % steer_vect: complex Nx1 (Number of elements)x1
            obj.received_signal = received_signal;  % Received signal at the sensor array
            obj.tx_num = tx_num;  % Number of transmitters
            obj.lambda = lambda;  % Wavelength (m)
            obj.element_num = element_num;  % Number of elements in the ULA
            obj.element_spacing = element_spacing;  % Distance between antenna elements (normalised to lambda units))
            obj.sweeping_angle = sweeping_angle;  % Angle range for sweeping to find the AoA
            obj.steer_vect = zeros(element_num, length(obj.sweeping_angle));
            obj.act_aoa = act_aoa;
            for i = 1:length(obj.sweeping_angle)
                obj.steer_vect(:, i) = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.sweeping_angle(i)) / obj.lambda);
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

        function square_err = calculate_square_error(obj, est_aoa)
            % Calculate the square error between the estimated and true AoA
            square_err = (est_aoa - obj.act_aoa).^2;
        end

        function result = parse_output(obj, type, spectrum)
            % Convert the spectrum to dB scale
            result.spectrum_dB = 10 * log10(abs(spectrum));
            % Find the peaks in the MVDR spectrum
            % [~, peak_indices] = findpeaks(spectrum_dB,'SortStr','descend');
            % est_aoa = obj.sweeping_angle(peak_indices(1:obj.tx_num));
            [~, max_idx] = max(result.spectrum_dB);
            result.est_aoa = obj.sweeping_angle(max_idx);
            result.square_err = obj.calculate_square_error(result.est_aoa);
            % obj.logger(type, est_aoa);
        end


        function result = ML_sync(obj, s_t)
            % Assume received signal: synchoronous signal with noise
            % Model: y = steering_vector*s_t + noise;
            % (know the transmitted symbols before channel and array effects)
            % s_t: complex 1xT (1 TX)x(Number of samples)
            % obj.received_signal: complex NxT (Number of elements)x(Number of samples)
            % steering_vector: complex Nx1 (Number of elements)x1
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                for t = 1:size(s_t,2) % consider each time instance separately
                    spectrum(i) = spectrum(i)+ real(obj.received_signal(:,t)' * obj.steer_vect(:, i) * s_t(:,t));
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
            for i = 1:length(obj.sweeping_angle)
                for t = 1:size(s_t,2) % consider each time instance separately
                    spectrum(i) = spectrum(i) + obj.received_signal(:,t)' * obj.steer_vect(:, i) * obj.steer_vect(:, i)' * obj.received_signal(:,t);
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
            for i = 1:length(obj.sweeping_angle)
                % Form MVDR denominator matrix from noise subspace eigenvectors
                spectrum(i) = obj.steer_vect(:, i)' * R * obj.steer_vect(:, i);  %steering_vector' * inv(R) * steering_vector;
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
            for i = 1:length(obj.sweeping_angle)
                % Form MVDR denominator matrix from noise subspace eigenvectors
                denom = (obj.steer_vect(:, i)' / R) * obj.steer_vect(:, i); % steering_vector' * inv(R) * steering_vector;
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
            for i = 1:length(obj.sweeping_angle)
                % Form MUSIC denominator matrix from noise subspace eigenvectors
                % another implementation for the denominator:
                %(steering_vector' * (noise_subspace * noise_subspace') * steering_vector +eps(1));
                denom = sum(abs(noise_subspace' * obj.steer_vect(:, i)).^2, 1);
                denom = denom+eps(1); % add a small positive constant to prevent division by zero. 9.44 in [1]
                spectrum(i) = 1 ./ denom;
            end
            % Parse the output to readable format
            result = obj.parse_output('MUSIC', spectrum);
        end
    end
end
