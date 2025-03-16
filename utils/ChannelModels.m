classdef ChannelModels < handle
    properties
    end
    
    methods
        function obj = ChannelModels()
        end

%% ---------------------- Generic Utils -----------------------------------
    function logger(~, rx_pos, pos_tx, angles, purpose)
        disp('--------------------------------------------------------');
        disp(purpose);
        disp('--------------------------------------------------------');
        disp(array2table(...
            angles, ...% table data
            'RowNames', cellstr(strcat('RX', num2str((1:size(rx_pos, 1))'))), ...
            'VariableNames', cellstr(strcat('TX', num2str((1:size(pos_tx, 1))')))));
    end

%% ---------------------- Position Characteristics ------------------------
    function [aoa_act, dist_act] = calculate_true_AoA_and_dist(~, rx_pos, pos_tx)
        aoa_act = zeros(size(rx_pos, 1), size(pos_tx, 1));
        dist_act = zeros(size(rx_pos, 1), size(pos_tx, 1));
        for i = 1:size(rx_pos, 1)
            for j = 1:size(pos_tx, 1)
                aoa_act(i,j) = atan2d(pos_tx(j,2) - rx_pos(i,2), pos_tx(j,1) - rx_pos(i,1)); % AoA in degrees - atan(y_tx-y_rx/x_tx-x_rx)
                dist_act = sqrt((pos_tx(j,1) - rx_pos(i,1))^2 + (pos_tx(j,2) - rx_pos(i,2))^2); % Euclidean distance between Tx and Rx - sqrt((x_tx-x_rx)^2 + (y_tx-y_rx)^2)
            end
        end
        % obj.logger(rx_pos, pos_tx, [aoa_act, dist_act], 'True Angles of Arrival with Absolute Distance');
    end
%% ---------------------- Channel Characteristics -------------------------
        function corrupted_output = AWGN(~, input, sigma_n2)
            % AWGN (independent among antenna elements)
            % Number of elements (N)
            % Number of samples (T)
            % Noise Variance (sigma_n2)
            [N, T] = size(input);
            Noise = sqrt(sigma_n2/2) * (randn(N, T) + 1j * randn(N, T));
            corrupted_output = input + Noise;
        end
        function output = FriisModel(~, input, distance, lambda)
            % Friis free space model which is dependent on distance
            % input: complex NxT (Number of TXs)x(Number of samples)
            % output: complex NxT (Number of TXs)x(Number of samples)
            % --- Channel gain using Friis free space model
            alpha = (lambda ./ (4 .* pi .* distance));
            % --- Apply channel gain to the transmitted signal
            output = alpha * input;
        end
            
        function output = LoS(~, input, amp_gain)
            % Light of Sight (LoS) model with a specified channel gain
            % input: complex NxT (Number of TXs)x(Number of samples)
            % output: complex NxT (Number of TXs)x(Number of samples)
            output = amp_gain .* input;
        end

        
        function gamma = computeGain(~, x_tx, y_tx, x_rx, y_rx, L_d0, d0, alpha)
            % Computes the steering vector gamma at each grid point.
            gamma = L_d0.^(-1/2) .* d0.^(alpha/2) .* ((x_rx-x_tx)^2+(y_rx-y_tx)^2).^(-alpha/4);
        end

%% ---------------------- Antenna Array Characteristics -------------------
        function output = applyULA(~, input, aoa_act, element_num, element_spacing, lambda)
            % ULA - Uniform Linear Array, we can apply the
            % input: complex 1xT (1 TXs)x(Number of samples)
            % output: complex NxT (Number of elements)x(Number of samples)
            % Normalised lambda to 1 for simplicity (i.e., setting fc=c=1)
            % --- steering vector to reflect the array characteristics
            % original steering vector: exp(-1j * 2 * pi * element_spacing * (0:N-1)' * sind(theta) / lambda)
            Alpha = exp(-1j * 2 * pi * element_spacing * (0:element_num-1)' * sind(aoa_act) / lambda);
            output = Alpha * input;
        end

%% ---------------------- Performance Evaluation --------------------------
        function CRB = CRB_det_1d(~, tx_sig, nPower, aoa_act, element_num, element_spacing, lambda)
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
            A = exp(-2j * pi * element_spacing * (0:element_num-1)' * sind(aoa_act) / lambda);
            D = A .* (-2j * pi * element_spacing * (0:element_num-1)' * cosd(aoa_act) / lambda);
            % --- Calculate the covariance matrix - sample correlation matrix
            P_est = tx_sig * tx_sig' / snapshot_count; % y*y^H/N
            [m, k] = size(A);
            H = D'*(eye(m) - A/(A'*A)*A')*D;
            CRB = real(H .* P_est.');
            CRB = eye(k) / CRB * (nPower / snapshot_count / 2);
        end

        function CRB = CRB_det_1d_simp(~, tx_sig, nPower, aoa_act, element_num, lambda)
            % --- CRB for a deterministic (conditional) model, in degree.
            % with simpler assumptions for the signal model.
            %
            % Inputs:
            %   tx_sig - Transmitted signal.
            %   nPower - Noise power.
            % There is no protection method yet
            % Implementation:
            % A = exp(-2j * pi * element_spacing * (0:element_num-1)' * sind(aoa_act) / lambda);
            % D = A .* (-2j * pi * element_spacing * (0:element_num-1)' * cosd(aoa_act) / lambda);
            % D2 = A.*(2j * pi * element_spacing * (0:element_num-1)' * sind(aoa_act) / lambda) + D .* (-2j * pi * element_spacing * (0:element_num-1)' * cosd(aoa_act) / lambda); 
            % % --- Calculate the covariance matrix - sample correlation matrix
            % CRB = -nPower / (tx_sig * tx_sig') / real(A' * D2);
            % Alternatively:
            CRB = 3/2 * nPower * lambda^2 / (tx_sig * tx_sig') / (element_num-1) / element_num / (2*element_num-1) / (pi*lambda/2*cosd(aoa_act)).^2;
        end

%% ---------------------- Generate received signal --------------------------
        function [nPower, y_centralised] = generateReceivedSignal(obj, sig_tx, pos_tx, pos_rx, aoa_act, e_avg, snr_db, L_d0, d0, alpha, element_num, element_spacing, lambda)
            % --- Generate received signal
            % Inputs:
            %   pos_tx - Position of the TXs.
            %   pos_rx - Position of the RXs.
            %   sig_tx - Transmitted signal.
            %   e_avg - Average energy of the transmitted signal.
            %   snr_db - SNR in dB.
            %   L_d0 - Path loss at reference distance d0.
            %   d0 - Reference distance.
            %   alpha - Path loss exponent.
            %   element_num - Number of elements in the array.
            %   element_spacing - Spacing between the elements.
            %   lambda - Wavelength.
            % Outputs:
            %   nPower - Noise power.
            %   y_centralised - Received signal.
            nvar_snr = length(snr_db);  % Number of SNR variants to test
            num_rx = size(pos_rx, 1);    % Update num_rx based on number of receivers
            y_centralised = cell(nvar_snr, num_rx); % Received signal at each Rx vectorised to cell array
            for idx_snr=1:nvar_snr
                % --- Generate received signal at each Rx
                for idx_rx=1:num_rx
                    nPower = e_avg/db2pow(snr_db(idx_snr, idx_rx));
                    y_los = obj.LoS(sig_tx, obj.computeGain(pos_tx(1), pos_tx(2), pos_rx(idx_rx, 1), pos_rx(idx_rx, 2), L_d0, d0, alpha));
                    y_ula = obj.applyULA(y_los, aoa_act(idx_rx), element_num, element_spacing, lambda);
                    y_awgn = obj.AWGN(y_ula, nPower);
                    y_centralised{idx_snr, idx_rx} = y_awgn;
                end
            end
        end
    end
end

