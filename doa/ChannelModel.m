classdef ChannelModel < handle
    properties
        tx_pos      % Transmitter position (x, y) in meters
        rx_pos      % Receiver position (x, y) in meters
        element_num % Number of elements in the ULA
        element_spacing    % Distance between antenna elements (normalised to lambda units)
        element_pos % Element relative positions (normalised to 0 at the first element)
        lambda      % Wavelength (m)
        aoa_act     % Actual Angle of Arrival
        dist_act    % Actual Distance
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
            [obj.aoa_act, obj.dist_act] = obj.calculate_true_AoA_and_dist();
        end
        
        function [aoa_act, dist_act] = calculate_true_AoA_and_dist(obj)
            aoa_act = zeros(size(obj.rx_pos, 1), size(obj.tx_pos, 1));
            dist_act = zeros(size(obj.rx_pos, 1), size(obj.tx_pos, 1));
            for i = 1:size(obj.rx_pos, 1)
                for j = 1:size(obj.tx_pos, 1)
                    aoa_act(i,j) = atan2d(obj.tx_pos(j,2) - obj.rx_pos(i,2), obj.tx_pos(j,1) - obj.rx_pos(i,1)); % AoA in degrees - atan(y_tx-y_rx/x_tx-x_rx)
                    dist_act = sqrt((obj.tx_pos(j,1) - obj.rx_pos(i,1))^2 + (obj.tx_pos(j,2) - obj.rx_pos(i,2))^2); % Euclidean distance between Tx and Rx - sqrt((x_tx-x_rx)^2 + (y_tx-y_rx)^2)
                end
            end
            % obj.logger('True Angles of Arrival with Absolute Distance', [aoa_act, dist_act]);
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
            % --- steering vector to reflect the array characteristics
            % original steering vector: exp(-1j * 2 * pi * element_spacing * (0:N-1)' * sind(theta) / lambda)
            Alpha = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.aoa_act) / obj.lambda);
            output = Alpha * input;
        end

        %% --- Performance Evaluation
        function CRB = CRB_det_1d(obj, tx_sig, nPower)
            % --- CRB for general 1D arrays based on the deterministic (conditional) model, in degree.
            %
            % Inputs:
            %   tx_sig - Transmitted signal.
            %   nPower - Noise power.
            %Reference:
            %   [1] P. Stoica and A. Nehorai, "Performance study of conditional and
            %       unconditional direction-of-arrival estimation," IEEE Transactions
            %       on Acoustics, Speech and Signal Processing, vol. 38, no. 10,
            %       pp. 1783-1795, Oct. 1990.
            snapshot_count = size(tx_sig, 2);
            A = exp(-2j * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.aoa_act) / obj.lambda);
            D = A .* (-2j * pi * obj.element_spacing * (0:obj.element_num-1)' * cosd(obj.aoa_act) / obj.lambda);
            % --- Calculate the covariance matrix - sample correlation matrix
            P_est = tx_sig * tx_sig' / snapshot_count; % y*y^H/N
            [m, k] = size(A);
            H = D'*(eye(m) - A/(A'*A)*A')*D;
            CRB = real(H .* P_est.');
            CRB = eye(k) / CRB * (nPower / snapshot_count / 2);
        end

        function CRB = CRB_det_1d_simp(obj, tx_sig, nPower)
            % --- CRB for a deterministic (conditional) model, in degree.
            % with simpler assumptions for the signal model.
            %
            % Inputs:
            %   tx_sig - Transmitted signal.
            %   nPower - Noise power.
            % There is no protection method yet
            % Implementation:
            % A = exp(-2j * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.aoa_act) / obj.lambda);
            % D = A .* (-2j * pi * obj.element_spacing * (0:obj.element_num-1)' * cosd(obj.aoa_act) / obj.lambda);
            % D2 = A.*(2j * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.aoa_act) / obj.lambda) + D .* (-2j * pi * obj.element_spacing * (0:obj.element_num-1)' * cosd(obj.aoa_act) / obj.lambda); 
            % % --- Calculate the covariance matrix - sample correlation matrix
            % CRB = -nPower / (tx_sig * tx_sig') / real(A' * D2);
            % Alternatively:
            CRB = 3/2 * nPower * obj.lambda^2 / (tx_sig * tx_sig') / (obj.element_num-1) / obj.element_num / (2*obj.element_num-1) / (pi*obj.lambda/2*cosd(obj.aoa_act)).^2;
        end
    end
end