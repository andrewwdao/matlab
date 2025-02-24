classdef Map3D < handle
    properties
        X  % Grid X values
        Y  % Grid Y values
        L  % Maximum Likelihood function values
    end
    
    methods
        function obj = Map3D(varargin)
            % Constructor: can take 1 or 3 arguments.
            % If 1 argument is provided, it's the ML matrix L and default X and Y grids are generated.
            % If 3 arguments are provided, they are X, Y, and L respectively.
            if nargin == 1
                obj.L = varargin{1};
                x = linspace(0, 100, 100);
                y = linspace(0, 100, 100);
                [obj.X, obj.Y] = meshgrid(x, y);
            elseif nargin == 3
                obj.X = varargin{1};
                obj.Y = varargin{2};
                obj.L = varargin{3};
            else
                error('Invalid number of arguments. Provide either 1 argument (likelihood matrix) or 3 arguments (X, Y, and likelihood matrix).');
            end
        end
        
        function plot(obj, ax)
            % Plots the 3D ML function on the provided axes.
            if nargin < 2 || isempty(ax)
                figure;
                ax = axes;
            end
            surf(ax, obj.X, obj.Y, abs(obj.L));
            xlabel(ax, 'x (m)');
            ylabel(ax, 'y (m)');
            zlabel(ax, 'L(x,y)');
            title(ax, '3D Visualization of L(x,y) = P(\theta_1, \theta_2)');
            shading(ax, 'interp'); % Improves plot appearance
            colorbar(ax); % Shows color scale
            view(ax, 30, 30);
        end
        
        function updateVisualization(obj, newX, newY, newL)
            % Update the grid and ML values then replot.
            obj.X = newX;
            obj.Y = newY;
            obj.L = newL;
        end
    end
end