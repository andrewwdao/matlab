classdef Channel
    properties
        lambda      % Wavelength (m), assume all Tx's operating at one centre freq at the moment
        distance    % Distances between Tx's and Rx (m) as a column vector, assume one Rx
        sigma_n2    % Noise variance (white noise)
        num_elements    % Number of antenna elements
        aoa      % Angles of arrival as a column vector, should have the same dimensionality as distance
        spacing_ula    % Antenna element spacing
    end
    
    methods
        function obj = Channel(lambda, distance, sigma_n2, num_elements, aoa, spacing_ula)
            % Constructor to initialize the properties
            obj.lambda = lambda;
            obj.distance = distance;
            obj.sigma_n2 = sigma_n2;
            obj.num_elements = num_elements;   
            obj.aoa = aoa; 
            obj.spacing_ula = spacing_ula; 
        end
        
       function Signal_before_antenna_array = LoS_specified_Loss(obj, transmitted_signal, amplitude_gain)
            % LoS with a specified channel gain
            % transmitted_signal: num_txs * num_samples complex signal
            % Signal_before_antenna_array: num_txs * num_samples complex signal with channel effects
            
            % Apply channel gain to the transmitted signal
            Signal_before_antenna_array = amplitude_gain .* transmitted_signal;    % dimension should be num_txs * num_timesamples
        end
        
        function Signal_before_antenna_array = FriisModel(obj, transmitted_signal)
            % Apply the Friis free space model
            % transmitted_signal: num_txs * num_samples complex signal
            % Signal_before_antenna_array: num_txs * num_samples complex signal with channel effects
            
            % Channel gain using Friis free space model
            friis_amplitude_gain = (obj.lambda ./ (4 * pi * obj.distance)) .* exp(-1j * 2 * pi * obj.distance / obj.lambda);
            
            % Apply channel gain to the transmitted signal
            Signal_before_antenna_array = friis_amplitude_gain .* transmitted_signal;    % dimension should be num_txs * num_timesamples
        end
        
        function Signal_after_antenna_array = ULA(obj, Signal_before_antenna_array)
            % ULA
            % steering matrix
            Alpha = exp(-1j * (0:(obj.num_elements-1))' * 2*pi * obj.spacing_ula * sind(obj.aoa)' / obj.lambda);   % dimension should be num_elements * num_txs
            % Multiply the rcvd signal by the rx steering matrix
            Signal_after_antenna_array = Alpha * Signal_before_antenna_array;   % dimension should be num_elements * num_timesamples
        end
        
        function AWGN_corrupted_Signal = AWGN(obj, Signal_after_antenna_array)
            % AWGN, also independent among antenna elements
            [N, T] = size(Signal_after_antenna_array);  % Number of elements (N) and number of samples (T)
            Noise = sqrt(obj.sigma_n2/2) * (randn(N, T) + 1j * randn(N, T));
            AWGN_corrupted_Signal = Signal_after_antenna_array + Noise;
        end
        
        
    end
end