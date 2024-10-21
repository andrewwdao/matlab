classdef ChannelModel < handle
    properties
        tx_pos      % Transmitter position (x, y) in meters
        rx_pos      % Receiver position (x, y) in meters
        element_num % Number of elements in the ULA
        element_spacing    % Distance between antenna elements (normalised to lambda units)
        element_pos % Element relative positions (normalised to 0 at the first element)
        lambda      % Wavelength (m)
    end
    
    methods
%% ------------------ Class initialisation and helpers -------------------
        function obj = ChannelModel(tx_pos, rx_pos, lambda, element_num, element_spacing)
            % Constructor
            obj.tx_pos = tx_pos;
            obj.rx_pos = rx_pos;
            obj.lambda = lambda;  % Wavelength (m)
            obj.element_num = element_num;  % Number of elements in the ULA
            obj.element_spacing = element_spacing;  % Distance between antenna elements (normalised to lambda units))
            obj.element_pos = 0:element_spacing:(element_num-1);  % Element relative positions (normalised to 0 at the first element)
        end
        
        function [act_aoa, act_dist] = calculate_true_AoA_and_dist(obj)
            act_aoa = zeros(size(obj.rx_pos, 1), size(obj.tx_pos, 1));
            act_dist = zeros(size(obj.rx_pos, 1), size(obj.tx_pos, 1));
            for i = 1:size(obj.rx_pos, 1)
                for j = 1:size(obj.tx_pos, 1)
                    act_aoa(i,j) = atan2d(obj.tx_pos(j,2) - obj.rx_pos(i,2), obj.tx_pos(j,1) - obj.rx_pos(i,1)); % AoA in degrees - atan(y_tx-y_rx/x_tx-x_rx)
                    act_dist = sqrt((obj.tx_pos(j,1) - obj.rx_pos(i,1))^2 + (obj.tx_pos(j,2) - obj.rx_pos(i,2))^2); % Euclidean distance between Tx and Rx - sqrt((x_tx-x_rx)^2 + (y_tx-y_rx)^2)
                end
            end
            % obj.logger('True Angles of Arrival with Absolute Distance', [act_aoa, act_dist]);
        end

        function logger(obj, purpose, angles)
            disp('--------------------------------------------------------');
            disp(purpose);
            disp('--------------------------------------------------------');
            disp(array2table(...
                angles, ...% table data
                'RowNames', cellstr(strcat('RX', num2str((1:size(obj.rx_pos, 1))'))), ...
                'VariableNames', cellstr(strcat('TX', num2str((1:size(obj.tx_pos, 1))')))));
        end
 %% ---------------------- Channel Characteristics -----------------------       
        function corrupted_output = AWGN(~, input, sigma_n2)
            % AWGN (independent among antenna elements)
            % Number of elements (N)
            % Number of samples (T)
            % Noise Variance (sigma_n2)
            [N, T] = size(input);
            Noise = sqrt(sigma_n2/2) * (randn(N, T) + 1j * randn(N, T));
            corrupted_output = input + Noise;
        end
        function output = FriisModel(obj, input, distance)
            % Friis free space model which is dependent on distance
            % input: complex NxT (Number of TXs)x(Number of samples)
            % output: complex NxT (Number of TXs)x(Number of samples)
            % --- Channel gain using Friis free space model
            alpha = (obj.lambda ./ (4 .* pi .* distance));
            % --- Apply channel gain to the transmitted signal
            output = alpha * input;
        end
            
        function output = LoS(~, input, amp_gain)
            % Light of Sight (LoS) model with a specified channel gain
            % input: complex NxT (Number of TXs)x(Number of samples)
            % output: complex NxT (Number of TXs)x(Number of samples)
            output = amp_gain .* input;
        end

%% -------------------- Antenna Array Characteristics --------------------
        function output = applyULA(obj, input)
            % ULA - Uniform Linear Array, we can apply the
            % input: complex 1xT (1 TXs)x(Number of samples)
            % output: complex NxT (Number of elements)x(Number of samples)
            % Normalised lambda to 1 for simplicity (i.e., setting fc=c=1)
            % --- calculating the actual angle of arrival from Tx and Rx locations
            [act_aoa, ~] = obj.calculate_true_AoA_and_dist();
            % --- steering vector to reflect the array characteristics
            % original steering vector: exp(-1j * 2 * pi * element_spacing * (0:N-1)' * sind(theta) / lambda)
            Alpha = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(act_aoa) / obj.lambda);
            output = Alpha * input';
        end
    end
end