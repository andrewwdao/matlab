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
            addParameter(p, 'aoa_est', {}, @iscell);
            addParameter(p, 'angle_limit', 60, @isnumeric);
            addParameter(p, 'show_options', [true, true]);
            addParameter(p, 'ax', [], @ishandle);
            addParameter(p, 'ShowAnnotation', false);
            addParameter(p, 'AnnotationPosition', [0.47, 0.81, 0.25, 0.10]);
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
            
            % Create 3D plot first to establish color limits
            subplot(1,2,1);
            surf(x_data, y_data, z_data, 'EdgeColor', 'none', 'FaceColor', 'interp');
            xlabel(p.Results.XLabel);
            ylabel(p.Results.YLabel);
            zlabel(p.Results.ZLabel);
            title('Likelihood Function');
            shading interp;
            colorbar;
            view(0, 90); % Top-down view
            
            % Store color limits for consistent scale - using modern 'clim' instead of 'caxis'
            clim_vals = clim;
            
            if show_map_view
                % Right Subplot: Map View with contours
                subplot(1,2,2); hold on;
                % Plot the 2D map first
                obj.map2d.plots(pos_tx, pos_rx, rot_abs, area_size, aoa_act, angle_limit, show_options);
                hold on;
                % Overlay contours from the likelihood function
                contour(x_data, y_data, z_data, 15, 'LineWidth', 0.8);

                % The rest of your contour handling can stay mostly the same
                % Send contours to back layer so they don't hide map elements
                contour_patches = findobj(gca, 'Type', 'patch');
                for i = 1:length(contour_patches)
                    % Set Z-data to send to back (below other elements)
                    z_data = get(contour_patches(i), 'ZData');
                    set(contour_patches(i), 'ZData', ones(size(z_data))*-1);
                end

                % Also handle the contour lines
                contour_lines = findobj(gca, 'Type', 'contour');
                if ~isempty(contour_lines)
                    uistack(contour_lines, 'bottom');
                end

                % Match color scaling with 3D plot
                clim(clim_vals);
                
                % Bring important elements to front
                uistack(findobj(gca, 'Type', 'line'), 'top');
                uistack(findobj(gca, 'Type', 'scatter'), 'top');
                uistack(findobj(gca, 'Type', 'text'), 'top');
            end
            
            % Add annotation if needed
            if ~isempty(AnnotationStrings)
                annotation('textbox', annotPos, ...
                    'String', AnnotationStrings, ...
                    'FitBoxToText', 'on', ...
                    'BackgroundColor', 'white', ...
                    'EdgeColor', 'black');
            end
        end
    end
end