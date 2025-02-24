classdef ULA < AntennaArray
    properties
        lambda          % Wavelength (m)
        element_num     % Number of elements in the ULA
        element_spacing % Distance between elements (normalized to lambda)
    end
    
    methods
        function obj = ULA(lambda, element_num, element_spacing)
            % Constructor: Initialize ULA parameters
            obj.lambda = lambda;
            obj.element_num = element_num;
            obj.element_spacing = element_spacing;
        end
        
        function steer_vect = getSteeringVector(obj, angles)
            % Compute steering vectors for given angles
            % Inputs:
            %   angles: Scalar or Nx1 vector of theta angles (degrees)
            % Outputs:
            %   steer_vect: element_num x N matrix of steering vectors
            if isrow(angles)
                angles = angles'; % Convert row vector to column vector
            end
            N = size(angles, 1);
            steer_vect = zeros(obj.element_num, N);
            for i = 1:N
                theta = angles(i, 1);
                phase = -2j * pi * obj.element_spacing * (0:obj.element_num-1)' * sind(theta) / obj.lambda;
                steer_vect(:, i) = exp(phase);
            end
        end
    end
end
