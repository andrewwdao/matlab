classdef AntennaArray < handle
    % Abstract base class for antenna array configurations
    methods (Abstract)
        steer_vect = getSteeringVector(obj, angles)
        % Inputs:
        %   angles: Matrix where each row is a set of angles
        %           (e.g., [theta] for ULA, [azimuth, elevation] for UCA/URA)
        % Outputs:
        %   steer_vect: Matrix of steering vectors (num_elements x num_angles)
    end
end
