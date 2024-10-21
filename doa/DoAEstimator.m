classdef DoAEstimator < handle
    properties
        received_signal % Received signal at the sensor array
        tx_num % Number of transmitters
        lambda % Wavelength (m)
        element_num % Number of elements in the ULA
        element_spacing % Distance between antenna elements (normalised to lambda units)
        sweeping_angle % Angle range for sweeping to find the AoA
    end

    methods
        function obj = DoAEstimator(received_signal, tx_num, lambda, element_num, element_spacing, sweeping_angle)
            % --- Constructor to initialize the properties
            obj.received_signal = received_signal;  % Received signal at the sensor array
            obj.tx_num = tx_num;  % Number of transmitters
            obj.lambda = lambda;  % Wavelength (m)
            obj.element_num = element_num;  % Number of elements in the ULA
            obj.element_spacing = element_spacing;  % Distance between antenna elements (normalised to lambda units))
            obj.sweeping_angle = sweeping_angle;  % Angle range for sweeping to find the AoA
        end

        function [est_aoa, spectrum_dB] = parse_output(obj, purpose, spectrum)
            % Convert the spectrum to dB scale
            spectrum_dB = 10 * log10(abs(spectrum));
            % Find the peaks in the MVDR spectrum
            % [~, peak_indices] = findpeaks(spectrum_dB,'SortStr','descend');
            % est_aoa = obj.sweeping_angle(peak_indices(1:obj.tx_num));
            [~, max_idx] = max(spectrum_dB);
            est_aoa = obj.sweeping_angle(max_idx);
            % obj.logger(purpose, est_aoa);
        end

        function [est_aoa, spectrum_dB] = ML_sync(obj)
            % Assume received signal: synchoronous signal with noise
            % y = steering_vector + noise;
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                steering_vector = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.sweeping_angle(i)) / obj.lambda);
                spectrum(i) = real(steering_vector' * obj.received_signal);
            end
            % Parse the output to readable format
            [est_aoa, spectrum_dB] = obj.parse_output('Simple Sync Maximum Likelihood', spectrum);
        end

        function [est_aoa, spectrum_dB] = ML_async(obj)
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                steering_vector = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.sweeping_angle(i)) / obj.lambda);
                % Form MVDR denominator matrix from noise subspace eigenvectors
                spectrum(i) = abs(steering_vector' * obj.received_signal);
            end
            % Parse the output to readable format
            [est_aoa, spectrum_dB] = obj.parse_output('Simple Sync Maximum Likelihood', spectrum);
        end

        function [est_aoa, spectrum_dB] = BF(obj)
            % --- Convensional Beamforming Algorithm
            % Similar to the ML estimator for the assumed asynchoronous received signal:
            % y = steering_vector*exp(j\theta) + noise;
            % --- Calculate the covariance matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2); % y*y^H/N
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                steering_vector = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.sweeping_angle(i)) / obj.lambda);
                % Form MVDR denominator matrix from noise subspace eigenvectors
                spectrum(i) = steering_vector' * R * steering_vector;  %steering_vector' * inv(R) * steering_vector;
            end
            % Parse the output to readable format
            [est_aoa, spectrum_dB] = obj.parse_output('Conventional Beamforming', spectrum);
        end
        
        function [est_aoa, spectrum_dB] = MVDR(obj)
            % --- MVDR Algorithm - or Minimum Variance Distortionless Response - Maximum Likelihood
            % --- Calculate the covariance matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2); % y*y^H/N
            % Add regularization to the covariance matrix (diagonal loading)
            R = R + 1e-3 * eye(size(R));  % add a small positive constant to prevent division by zero. 9.44 in [1]
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                steering_vector = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.sweeping_angle(i)) / obj.lambda);
                % Form MVDR denominator matrix from noise subspace eigenvectors
                denom = (steering_vector' / R) * steering_vector; % steering_vector' * inv(R) * steering_vector;
                spectrum(i) = 1 ./ denom;
            end
            % Parse the output to readable format
            [est_aoa, spectrum_dB] = obj.parse_output('MVDR', spectrum);
        end

        function [est_aoa2, spectrum_dB2] = MUSIC(obj)
            % --- MUSIC Algorithm
            % --- Calculate the covariance matrix - sample correlation matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2); % y*y^H/N
            % Perform eigenvalue decomposition
            [eigenvectors, eigenvalues] = eig(R);
            % Sort eigenvalues and eigenvectors
            [~, idx] = sort(diag(eigenvalues), 'descend');
            eigenvectors = eigenvectors(:, idx); % get the largest eigenvectors
            % Determine the noise subspace
            noise_subspace = eigenvectors(:, obj.tx_num+1:end);
            % Compute the MUSIC spectrum
            spectrum = zeros(size(obj.sweeping_angle));
            for i = 1:length(obj.sweeping_angle)
                steering_vector = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.sweeping_angle(i)) / obj.lambda);
                % Form MUSIC denominator matrix from noise subspace eigenvectors
                % another implementation for the denominator:
                %(steering_vector' * (noise_subspace * noise_subspace') * steering_vector +eps(1));
                denom = sum(abs(noise_subspace' * steering_vector).^2, 1)+eps(1); % add a small positive constant to prevent division by zero. 9.44 in [1]
                spectrum(i) = 1 ./ denom;
            end
            % Parse the output to readable format
            [est_aoa2, spectrum_dB2] = obj.parse_output('MUSIC', spectrum);
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
    end
end