classdef PosEstimator2D < handle
    properties
    end
    methods
        function obj = PosEstimator2D()
        end

        function abs_ray = calAbsRay(~, rx_abs_pos, tx_abs_pos, rx_abs_rot, rel_angle)
            % Calculate the absolute ray parameters (slope and shift)
            % for a line connecting the Tx and Rx in the 2D plane.
            %
            % This is an absolute slope and shift added to
            % a relative angle of arrival for an RX at the origin (0,0) with no shift.
            % 
            % The function also determines the limits of the ray based on the quadrant
            % in which the Rx and Tx are located.
            %
            % Inputs:
            %    rx_abs_pos - A 1x2 vector representing the absolute position of the Rx in the 2D plane [x, y].
            %    tx_abs_pos - A 1x2 vector representing the absolute position of the Tx in the 2D plane [x, y].
            %    rx_abs_rot - A scalar representing the absolute rotation of the Rx in degrees.
            %    rel_angle  - A scalar representing the relative angle of arrival from TX to RX in degrees
            %                 that is output from the DoA Algorithms.
            %
            % Outputs:
            %    abs_ray    - A structure containing the following fields:
            %                 slope - The absolute slope of the line connecting the Tx and Rx.
            %                 shift - The absolute shift of the line connecting the Tx and Rx.
            %                 lim   - A 1x2 vector representing the limits of the ray in the 2D plane.
            %
            % Example:
            %    rx_abs_pos = [100, 200];
            %    tx_abs_pos = [300, 400];
            %    rx_abs_rot = 45;
            %    rel_angle = 30;
            %    abs_ray = calAbsRay([], rx_abs_pos, tx_abs_pos, rx_abs_rot, rel_angle);
            %
            % Other m-files required: None
            % Subfunctions: None
            % MAT-files required: None
            abs_ray.slope = tand(rel_angle + rx_abs_rot);
            abs_ray.shift = rx_abs_pos(2) - abs_ray.slope * rx_abs_pos(1);
            % abs_ray.lim = sort(abs_ray.slope * [0, 1000] + rx_abs_pos(1), 'ascend'); % Absolute ray limit from a relative x limit from 0 to 100
            if (rx_abs_pos(1) < tx_abs_pos(1) && rx_abs_pos(2) < tx_abs_pos(2)) || (rx_abs_pos(1) > tx_abs_pos(1) && rx_abs_pos(2) < tx_abs_pos(2))
                lim_range = [0, 1000]; % Third quadrant or first quadrant
            else
                lim_range = [-1000, 0]; % Second quadrant or fourth quadrant
            end
            abs_ray.lim = sort(abs_ray.slope * lim_range + rx_abs_pos(1), "ascend"); % Absolute ray limit from a relative x limit from 0 to 100
        end

        function result = calIntersection(~, abs_ray1, abs_ray2)
            % Calculate the intersection point of two lines in the 2D plane.
            % 
            % Inputs:
            %    abs_ray1 - The first line in the 2D plane.
            %    abs_ray2 - The second line in the 2D plane.
            %
            % Outputs:
            %    result - The intersection point of the two lines.
            %
            % Other m-files required: None
            % Subfunctions: None
            % MAT-files required: None
            result.x = (abs_ray2.shift - abs_ray1.shift) / (abs_ray1.slope - abs_ray2.slope);
            result.y = abs_ray1.slope * result.x + abs_ray1.shift;
        end
    end
end
