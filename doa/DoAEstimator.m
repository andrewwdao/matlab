classdef DoAEstimator < handle
    properties
        received_signal % Received signal at the sensor array
        tx_num % Number of transmitters
        lambda % Wavelength (m)
        element_num % Number of elements in the ULA
        element_spacing % Distance between antenna elements (normalised to lambda units)
        sweeping_angle % Angle range for sweeping to find the AoA
        sigma_n2 % Noise variance (white noise)
    end

    methods
        function obj = DoAEstimator(received_signal, tx_num, lambda, element_num, element_spacing, sweeping_angle, sigma_n2)
            % --- Constructor to initialize the properties
            obj.received_signal = received_signal;  % Received signal at the sensor array
            obj.tx_num = tx_num;  % Number of transmitters
            obj.lambda = lambda;  % Wavelength (m)
            obj.element_num = element_num;  % Number of elements in the ULA
            obj.element_spacing = element_spacing;  % Distance between antenna elements (normalised to lambda units))
            obj.sweeping_angle = sweeping_angle;  % Angle range for sweeping to find the AoA
            obj.sigma_n2 = sigma_n2;  % Noise variance (white noise)
        end

        function [est_aoa, music_spectrum_dB] = MUSIC(obj)
            % --- MUSIC Algorithm
            % --- Calculate the covariance matrix
            R = obj.received_signal * obj.received_signal' / size(obj.received_signal, 2); % y*y^H/N
            % Perform eigenvalue decomposition
            [eigenvectors, eigenvalues] = eig(R);
            % Sort eigenvalues and eigenvectors
            [eigenvalues, idx] = sort(diag(eigenvalues), 'descend');
            eigenvectors = eigenvectors(:, idx); % get the largest eigenvectors
            % Determine the noise subspace
            num_signals = 1; % Number of signals (assuming 1 for simplicity)
            noise_subspace = eigenvectors(:, num_signals+1:end);
            % Compute the MUSIC spectrum
            angles = -90:1:90; % Angle range for MUSIC spectrum
            music_spectrum = zeros(size(angles));

            for i = 1:length(angles)
                steering_vector = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(angles(i)) / obj.lambda);
                % Form MUSIC denominator matrix from noise subspace eigenvectors
                denom = sum(abs(noise_subspace' * steering_vector).^2, 1)+eps(1); % add a small positive constant to prevent division by zero. 9.44 in [1]
                % another implementation for the denominator:
                %(steering_vector' * (noise_subspace * noise_subspace') * steering_vector +eps(1));
                music_spectrum(i) = 1 ./ denom;
            end

            % Convert MUSIC spectrum to dB scale
            music_spectrum_dB = 10 * log10(abs(music_spectrum));

            % Find the peaks in the MUSIC spectrum
            [~, peak_indices] = findpeaks(music_spectrum_dB,'SortStr','descend');
            est_aoa = obj.sweeping_angle(peak_indices(1));
            obj.log('Estimated Angles of Arrival', est_aoa);
        end
        
        function log(obj, purpose, values)
            disp('--------------------------------------------------------');
            disp(purpose);
            disp('--------------------------------------------------------');
            disp(array2table(...
                values, ...% table data
                'RowNames', cellstr(strcat('RX', num2str((1:1)'))), ...
                'VariableNames', cellstr(strcat('TX', num2str((1:obj.tx_num)')))));
        end
    end
end