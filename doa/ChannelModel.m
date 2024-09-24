classdef ChannelModel
    properties
        lambda      % Wavelength (m)
        distance    % Distance between Tx and Rx (m)
        sigma_n2    % Noise variance (white noise)
    end
    
    methods
        function obj = ChannelModel(lambda, distance, sigma_n2)
            % Constructor to initialize the properties
            obj.lambda = lambda;
            obj.distance = distance;
            obj.sigma_n2 = sigma_n2;
        end
        
        function received_signal = FriisModel(obj, transmitted_signal)
            % Apply the Friis free space model with Gaussian white noise
            % transmitted_signal: NxT complex signal
            % received_signal: NxT complex signal with channel effects
            
            % Channel gain using Friis free space model
            alpha = (obj.lambda / (4 * pi * obj.distance)) * exp(-1j * 2 * pi * obj.distance / obj.lambda);
            
            % Apply channel gain to the transmitted signal
            received_signal = alpha * transmitted_signal;
            
            % Add Gaussian white noise
            [N, T] = size(transmitted_signal);  % Number of elements (N) and number of samples (T)
            noise = sqrt(obj.sigma_n2) * (randn(N, T) + 1j * randn(N, T));
            received_signal = received_signal + noise;
        end
    end
end