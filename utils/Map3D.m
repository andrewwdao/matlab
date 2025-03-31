classdef Map3D < handle
    % Map3D class for visualizing 3D Maximum Likelihood (ML) function and 2D map view
    %   This class provides methods to plot the 3D ML function and optionally a 2D map view.
    %   It is designed to work with the Map2D class for 2D visualization.
    %
    properties
        map2d % Map2D object for 2D visualization
    end
    methods
        function obj = Map3D(map2d, varargin)
            %% Constructor for Map3D class
            obj.map2d = map2d; % Store the Map2D object
            % Parse additional parameters
        end
        
        function plots(obj, x_data, y_data, z_data, varargin)
            %% Plots the 3D ML function and optionally a 2D map view
            %
            % This function creates a 3D surface plot of likelihood data and can optionally display a 2D map view alongside it.
            % When a map view is shown, the figure is split into two subplots with the 3D visualization on the left
            % and the map view on the right.
            %
            % Usage:
            %   obj.plots(x_data, y_data, z_data)
            %   obj.plots(x_data, y_data, z_data, 'Name', Value, ...)
            %
            % Parameters:
            %   x_data          - X coordinates (meshgrid format)
            %   y_data          - Y coordinates (meshgrid format)
            %   z_data          - Likelihood values at each (x,y) coordinate
            %
            % Optional Name-Value Pairs:
            %   'Title'             - Plot title (default: '3D Visualization of The Likelihood Function')
            %   'XLabel'            - X-axis label (default: 'x (m)')
            %   'YLabel'            - Y-axis label (default: 'Error [m]')
            %   'ZLabel'            - Z-axis label (default: 'L(x,y)')
            %   'pos_tx'            - Transmitter position [x,y] (default: [])
            %   'pos_rx'            - Receiver positions matrix (default: [])
            %   'rot_abs'           - Absolute rotations (default: [])
            %   'area_size'         - Size of the area to plot (default: 100)
            %   'aoa_act'           - Actual angles of arrival (default: [])
            %   'angle_limit'       - Angle limit for visualization (default: 60)
            %   'show_options'      - Options for what to display [bool, bool] (default: [true, true])
            %   'ax'                - Axes handle to plot on (default: [])
            %   'AnnotationPosition'- Position of the annotation box [left, bottom, width, height] (default: [0.01, 0.85, 0.25, 0.10])
            %   'AnnotationStrings' - Cell array of strings for annotation (default: {})
            %
            % Note:
            %   The map view is only shown if obj.map2d is not empty and both pos_tx and pos_rx are provided.
            %   The 3D visualization is displayed with a top-down view (view(0, 90)).
            %   X, Y - Meshgrid coordinates
            %   L - Likelihood values
            %   Optional Name-Value Pairs:
            %     'AnnotationStrings' - Cell array of strings for annotation
            %     'pos_tx' - Transmitter position
            %     'pos_rx' - Receiver positions
            %     'rot_abs' - Absolute rotations
            %     'area_size' - Area size
            %     'aoa_act' - Actual angles of arrival
            %     'angle_limit' - Angle limit
            
            % Parse inputs
            p = inputParser;
            addParameter(p, 'Title', '3D Visualization', @ischar);
            addParameter(p, 'XLabel', 'x (m)', @ischar);
            addParameter(p, 'YLabel', 'y (m)', @ischar);
            addParameter(p, 'ZLabel', 'L(x,y)', @ischar);
            addParameter(p, 'pos_tx', [], @isnumeric);
            addParameter(p, 'pos_rx', [], @isnumeric);
            addParameter(p, 'rot_abs', [], @isnumeric);
            addParameter(p, 'area_size', 100, @isnumeric);
            addParameter(p, 'aoa_act', [], @isnumeric);
            addParameter(p, 'angle_limit', 60, @isnumeric);
            addParameter(p, 'show_options', [true, true]);
            addParameter(p, 'ax', [], @ishandle);
            addParameter(p, 'ShowAnnotation', false);
            addParameter(p, 'AnnotationPosition', [0.01, 0.85, 0.25, 0.10]);
            addParameter(p, 'AnnotationStrings', {}, @iscell);
            parse(p, varargin{:});
            
            annotPos = p.Results.AnnotationPosition;
            AnnotationStrings = p.Results.AnnotationStrings;
            pos_tx = p.Results.pos_tx;
            pos_rx = p.Results.pos_rx;
            rot_abs = p.Results.rot_abs;
            area_size = p.Results.area_size;
            aoa_act = p.Results.aoa_act;
            angle_limit = p.Results.angle_limit;
            show_options = p.Results.show_options;
            ax = p.Results.ax;
            
            % Create figure if needed
            if isempty(ax)
                figure('Name', p.Results.Title, 'WindowState', 'maximized', 'NumberTitle', 'off');
            end
            
            % Determine if we should show map view
            show_map_view = ~isempty(obj.map2d) && ~isempty(pos_tx) && ~isempty(pos_rx);
            
            if show_map_view
                % Right Subplot: Map View
                subplot(1,2,2); hold on;
                obj.map2d.plot(pos_tx, pos_rx, rot_abs, area_size, aoa_act, angle_limit, show_options);
                % Left Subplot: 3D Visualization of ML Function
                subplot(1,2,1);
            end
            
            % Create 3D plot
            surf(x_data, y_data, z_data, 'EdgeColor', 'none', 'FaceColor', 'interp');
            xlabel(p.Results.XLabel);
            ylabel(p.Results.YLabel);
            zlabel(p.Results.ZLabel);
            title(p.Results.Title);
            shading interp; % Improves plot appearance
            colorbar; % Shows color scale
            view(0, 90); % Top-down view
            
            % Add annotation if AnnotationStrings provided
            if p.Results.ShowAnnotation && ~isempty(AnnotationStrings)
                annotation('textbox', annotPos, ...
                    'String', AnnotationStrings, ...
                    'FitBoxToText', 'on', ...
                    'BackgroundColor', 'white', ...
                    'EdgeColor', 'black');
            end
        end
    end
end