classdef ArrayModel < handle
    properties
        tx_pos      % Transmitter position (x, y) in meters
        rx_pos      % Receiver position (x, y) in meters
        act_aoa    % True Angle of Arrival (AoA)
        act_distance % Actual distance between Tx and Rx
        element_num % Number of elements in the ULA
        element_spacing    % Distance between antenna elements (normalised to lambda units)
        element_pos % Element relative positions (normalised to 0 at the first element)
        lambda      % Wavelength (m)
    end
    
    methods
        function obj = ArrayModel(tx_pos, rx_pos, lambda, element_num, element_spacing)
            % --- Constructor to initialize the properties
            obj.element_num = element_num;  % Number of elements in the ULA
            obj.element_spacing = element_spacing;  % Distance between antenna elements (normalised to lambda units))
            obj.element_pos = 0:element_spacing:(element_num-1);  % Element relative positions (normalised to 0 at the first element)
            obj.tx_pos = tx_pos;
            % --- True Angle of Arrival (AoA)
            obj.act_aoa = zeros(size(rx_pos, 1), size(tx_pos, 1));
            for i = 1:size(rx_pos, 1)
                for j = 1:size(tx_pos, 1)
                    obj.act_aoa(i,j) = atan2d(tx_pos(j,2) - rx_pos(i,2), tx_pos(j,1) - rx_pos(i,1)); % AoA in degrees - atan(y_tx-y_rx/x_tx-x_rx)
                    obj.act_distance = sqrt((tx_pos(j,1) - rx_pos(i,1))^2 + (tx_pos(j,2) - rx_pos(i,2))^2); % Euclidean distance between Tx and Rx - sqrt((x_tx-x_rx)^2 + (y_tx-y_rx)^2)
                end
            end
            % obj.logger('True Angles of Arrival with Absolute Distance', sprintf("%.2f deg; %.2f m", obj.act_aoa, obj.act_distance));
            % obj.logger('True Angles of Arrival with Absolute Distance', obj.act_aoa);
            obj.lambda = lambda;  % Wavelength (m)
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

        function output = applyULA(obj, input)
            % Since this is a ULA - Uniform Linear Array, we can apply the
            % steering vector to the input signal to reflect the array
            % characteristics
            % input: NxT complex signal
            % output: NxT complex signal with array characteristics
            % Normalised lambda to 1 for simplicity (i.e., setting fc=c=1)
            % original steering vector: exp(-1j * 2 * pi * element_spacing * (0:N-1)' * sind(theta) / lambda)
            output = exp(-1j * 2 * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(obj.act_aoa) / obj.lambda) * input';
        end
    end
end