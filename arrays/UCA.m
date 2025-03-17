
% This is a simplified implementation.
% You can refine it based on your specific UCA model when needed.
classdef UCA < AntennaArray
    properties
        lambda
        radius          % Radius of the circle (in wavelengths)
        element_num
    end
    
    methods
        function obj = UCA(lambda, radius, element_num)
            obj.lambda = lambda;
            obj.radius = radius;
            obj.element_num = element_num;
        end
        
        function steer_vect = getSteeringVector(obj, angles)
            % Compute steering vectors for UCA
            % angles: Nx2 matrix of [azimuth, elevation] pairs (degrees)
            steer_vect = zeros(obj.element_num, size(angles, 1));
            phi = linspace(0, 2*pi, obj.element_num+1);
            phi = phi(1:end-1); % Angular positions of elements
            for i = 1:size(angles, 1)
                az = angles(i, 1); % Azimuth
                el = angles(i, 2); % Elevation
                phase = -2j * pi * (obj.radius / obj.lambda) * cos(phi - deg2rad(az)) .* sin(deg2rad(el));
                steer_vect(:, i) = exp(phase);
            end
        end
    end
end
