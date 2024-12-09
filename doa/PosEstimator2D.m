classdef PosEstimator2D < handle
    properties
    end
    methods
        function obj = PosEstimator2D()
        end

        function ray_plot_lim = calRayPlotLim(~, pos_abs_rx, pos_abs_tx, slope)
            % Calculate the limits of the ray plot based on the quadrant
            % in which the Rx and Tx are located.
            % 
            % Inputs:
            %    pos_abs_rx - A 1x2 vector representing the absolute position of the Rx in the 2D plane [x, y].
            %    pos_abs_tx - A 1x2 vector representing the absolute position of the Tx in the 2D plane [x, y].
            %    slope      - A scalar representing the slope of the line connecting the Tx and Rx.
            % 
            % Outputs:
            %    ray_plot_lim - A 1x2 vector representing the limits of the ray plot in the 2D plane.
            % 
            % Example:
            %    pos_abs_rx = [100, 200];
            %    pos_abs_tx = [300, 400];
            %    slope = 0.5;
            %    ray_plot_lim = calRayPlotLim(pos_abs_rx, pos_abs_tx, slope);
            % 
            % Other m-files required: None
            % Subfunctions: None
            % MAT-files required: None
            if (pos_abs_rx(1) < pos_abs_tx(1) && pos_abs_rx(2) < pos_abs_tx(2)) || (pos_abs_rx(1) > pos_abs_tx(1) && pos_abs_rx(2) < pos_abs_tx(2))
                lim_range = [0, 10000]; % Third quadrant or first quadrant
            else
                lim_range = [-10000, 0]; % Second quadrant or fourth quadrant
            end
            ray_plot_lim = sort(slope * lim_range + pos_abs_rx(1), "ascend"); % Absolute ray limit from a relative x limit from 0 to 100
        end

        function abs_ray = calAbsRays(obj, pos_abs_rx, pos_abs_tx, rot_abs_rx, aoa_rel, rot_lim)
            % Calculate the absolute ray parameters (slope and doa_shift)
            % for a line connecting the Tx and Rx in the 2D plane.
            %
            % This is an absolute slope and doa_shift added to
            % a relative angle of arrival for an RX at the origin (0,0) with no doa_shift.
            % 
            % The function also determines the limits of the ray based on the quadrant
            % in which the Rx and Tx are located.
            %
            % Inputs:
            %    pos_abs_rx - A 1x2 vector representing the absolute position of the Rx in the 2D plane [x, y].
            %    pos_abs_tx - A 1x2 vector representing the absolute position of the Tx in the 2D plane [x, y].
            %    rot_abs_rx - A scalar representing the absolute rotation of the Rx in degrees.
            %    aoa_rel - A scalar representing the relative angle of arrival from TX to RX in degrees
            %                 that is output from the DoA Algorithms.
            %    rot_lim (optional) - A scalar representing the limit of the rotation in degrees.
            %
            % Outputs:
            %    abs_ray    - A structure containing the following fields:
            %                 centre_slope - The slope of the line of the absolute rotation of the Rx. to the x-axis.
            %                 cw_slope     - The slope of the line of the clock-wise absolute rotation of the Rx to the x-axis.
            %                 ccw_slope    - The slope of the line of the counter-clock-wise absolute rotation of the Rx to the x-axis.
            %                 doa_slope    - The slope of the line of the absolute angle of arrival of the Rx to the Tx.
            %                 doa_shift - The absolute doa_shift of the lines that go through RX.
            %                 lim   - A 1x2 vector representing the limits of the ray in the 2D plane that is used for plotting.
            %
            % Example:
            %    pos_abs_rx = [100, 200];
            %    pos_abs_tx = [300, 400];
            %    rot_abs_rx = 45;
            %    aoa_rel = 30;
            %    abs_ray = calAbsRays([], pos_abs_rx, pos_abs_tx, rot_abs_rx, aoa_rel);
            %
            % Other m-files required: None
            % Subfunctions: None
            % MAT-files required: None
            abs_ray.doa_slope = tand(rot_abs_rx + aoa_rel);
            abs_ray.doa_shift = pos_abs_rx(2) - abs_ray.doa_slope * pos_abs_rx(1);
            if nargin == 6
                abs_ray.centre_slope = tand(rot_abs_rx);
                abs_ray.centre_shift = pos_abs_rx(2) - abs_ray.centre_slope * pos_abs_rx(1);
                abs_ray.cw_slope = tand(rot_abs_rx - rot_lim);
                abs_ray.cw_shift = pos_abs_rx(2) - abs_ray.cw_slope * pos_abs_rx(1);
                abs_ray.ccw_slope = tand(rot_abs_rx + rot_lim);
                abs_ray.ccw_shift = pos_abs_rx(2) - abs_ray.ccw_slope * pos_abs_rx(1);
                abs_ray.lim = obj.calRayPlotLim(pos_abs_rx, pos_abs_tx, abs_ray.doa_slope);
            end
        end

        function result = calDoAIntersect(~, abs_ray1, abs_ray2)
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
            result.x = (abs_ray2.doa_shift - abs_ray1.doa_shift) / (abs_ray1.doa_slope - abs_ray2.doa_slope);
            result.y = abs_ray1.doa_slope * result.x + abs_ray1.doa_shift;
        end
    end
end
