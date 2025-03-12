classdef Metric < handle
    % Metric - A class for calculating various error metrics in statistical analysis
    % 
    % This class provides methods for calculating common error metrics used in 
    % performance evaluation of estimation algorithms, such as MSE, RMSE, median 
    % error, and percentile-based metrics.
    
    properties
        percentile_values = [25, 50, 75];  % Default percentiles to calculate
    end
    
    methods
        function obj = Metric(varargin)
            % METRIC Constructor for the Metric class
            %   obj = METRIC() creates a Metric object with default percentile values
            %   obj = METRIC(percentile_values) creates a Metric object with specified 
            %   percentile values
            
            if nargin > 0
                obj.percentile_values = varargin{1};
            end
        end
        
        function mse = cal_MSE(~, errors)
            % CALCULATEMSE Calculate Mean Square Error
            %   mse = CALCULATEMSE(errors) calculates the mean square error from 
            %   the input error values
            %
            %   Input:
            %     errors: Cell array or matrix of error values
            %   Output:
            %     mse: Mean Square Error value(s)
            
            if iscell(errors)
                % For cell array of errors (e.g., from multiple SNRs/angles)
                [rows, cols] = size(errors);
                mse = zeros(rows, cols);
                
                for i = 1:rows
                    for j = 1:cols
                        mse(i, j) = mean(errors{i, j}.^2);
                    end
                end
            else
                % For direct array of errors
                mse = mean(errors.^2);
            end
        end
        
        function rmse = cal_RMSE(obj, errors)
            % CALCULATERMSE Calculate Root Mean Square Error
            %   rmse = CALCULATERMSE(errors) calculates the root mean square error 
            %   from the input error values
            %
            %   Input:
            %     errors: Cell array or matrix of error values
            %   Output:
            %     rmse: Root Mean Square Error value(s)
            
            mse = obj.cal_MSE(errors);
            rmse = sqrt(mse);
        end
        
        function percentiles = cal_Percentiles(~, errors, percentile_values)
            % CAL_PERCENTILES Calculate specified percentiles of error values
            %   percentiles = CAL_PERCENTILES(errors) calculates the 50th percentile (median)
            %   percentiles = CAL_PERCENTILES(errors, 50) also calculates the median
            %   percentiles = CAL_PERCENTILES(errors, [25, 50, 75]) calculates lower quartile,
            %   median, and upper quartile
            %
            %   Input:
            %     errors: Cell array or matrix of error values
            %     percentile_values: Single value or array of exactly 3 percentile values
            %   Output:
            %     percentiles: Structure with fields .val (for single percentile) or
            %                  .lower, .val, .upper (for 3 percentiles)
            %
            %   If percentile_values is not provided, defaults to median (50).
            %   Function will error if percentile_values contains anything other than
            %   1 or exactly 3 values.

            % Check and set default percentile values
            if nargin < 3 || isempty(percentile_values)
                percentile_values = 50;  % Default to median
            end

            % Validate percentile_values
            if ~(isscalar(percentile_values) || (isvector(percentile_values) && length(percentile_values) == 3))
                error('percentile_values must be either a single value or exactly 3 values');
            end

            % Initialize output struct based on number of percentiles
            is_single_percentile = isscalar(percentile_values);

            if iscell(errors)
                % For cell array of errors
                [rows, cols] = size(errors);
                
                % Initialize structure with appropriate fields
                if is_single_percentile
                    percentiles = struct('val', zeros(rows, cols));
                else
                    % For 3 percentiles, assume [lower, median, upper]
                    percentiles = struct('lower', zeros(rows, cols), ...
                                        'val', zeros(rows, cols), ...
                                        'upper', zeros(rows, cols));
                    % Sort percentile values to ensure lower/middle/upper order
                    percentile_values = sort(percentile_values);
                end
                
                % Calculate percentiles
                for i = 1:rows
                    for j = 1:cols
                        if is_single_percentile
                            percentiles.val(i, j) = prctile(errors{i, j}, percentile_values);
                        else
                            percentiles.lower(i, j) = prctile(errors{i, j}, percentile_values(1));
                            percentiles.val(i, j) = prctile(errors{i, j}, percentile_values(2));
                            percentiles.upper(i, j) = prctile(errors{i, j}, percentile_values(3));
                        end
                    end
                end
            else
                % For direct array of errors
                if is_single_percentile
                    percentiles = struct('val', prctile(errors, percentile_values));
                else
                    percentile_values = sort(percentile_values);
                    percentiles = struct('lower', prctile(errors, percentile_values(1)), ...
                                    'val', prctile(errors, percentile_values(2)), ...
                                    'upper', prctile(errors, percentile_values(3)));
                end
            end
        end
        
        function plotErrorMetrics(~, x_data, error_data, plot_type, varargin)
            % PLOTERRORMETRICS Create plot of error metrics
            %   PLOTERRORMETRICS(x_data, error_data, plot_type) creates a plot of error metrics
            %   with x_data on the x-axis and error_data on the y-axis
            %
            %   Input:
            %     x_data: Data for x-axis
            %     error_data: Data for y-axis (error metrics)
            %     plot_type: Type of plot ('linear', 'semilogy', 'loglog', etc.)
            %     varargin: Additional arguments for plot customization
            %
            %   Optional Name-Value Pairs:
            %     'LineStyle' - Line style for the plot
            %     'Marker' - Marker for the plot
            %     'LineWidth' - Width of the line
            %     'MarkerSize' - Size of the marker
            %     'Color' - Color of the line
            %     'DisplayName' - Name for the legend
            %     'ShowBand' - Whether to show the error band
            %     'BandLower' - Lower values for error band
            %     'BandUpper' - Upper values for error band
            
            % Parse inputs
            p = inputParser;
            addParameter(p, 'LineStyle', '-');
            addParameter(p, 'Marker', 'o');
            addParameter(p, 'LineWidth', 2);
            addParameter(p, 'MarkerSize', 6);
            addParameter(p, 'Color', []);
            addParameter(p, 'DisplayName', '');
            addParameter(p, 'ShowBand', false);
            addParameter(p, 'BandLower', []);
            addParameter(p, 'BandUpper', []);
            parse(p, varargin{:});
            
            % Create plot based on type
            switch lower(plot_type)
                case 'linear'
                    h_line = plot(x_data, error_data, ...
                        'LineStyle', p.Results.LineStyle, ...
                        'Marker', p.Results.Marker, ...
                        'LineWidth', p.Results.LineWidth, ...
                        'MarkerSize', p.Results.MarkerSize, ...
                        'DisplayName', p.Results.DisplayName);
                case 'semilogy'
                    h_line = semilogy(x_data, error_data, ...
                        'LineStyle', p.Results.LineStyle, ...
                        'Marker', p.Results.Marker, ...
                        'LineWidth', p.Results.LineWidth, ...
                        'MarkerSize', p.Results.MarkerSize, ...
                        'DisplayName', p.Results.DisplayName);
                case 'loglog'
                    h_line = loglog(x_data, error_data, ...
                        'LineStyle', p.Results.LineStyle, ...
                        'Marker', p.Results.Marker, ...
                        'LineWidth', p.Results.LineWidth, ...
                        'MarkerSize', p.Results.MarkerSize, ...
                        'DisplayName', p.Results.DisplayName);
                otherwise
                    error('Unknown plot type: %s', plot_type);
            end
            
            % Set color if specified
            if ~isempty(p.Results.Color)
                set(h_line, 'Color', p.Results.Color);
            end
            
            % Add confidence band if requested
            if p.Results.ShowBand && ~isempty(p.Results.BandLower) && ~isempty(p.Results.BandUpper)
                % Get line color
                line_color = get(h_line, 'Color');
                
                % Create shaded area
                hold on;
                fill([x_data; flipud(x_data)], ...
                    [p.Results.BandLower; flipud(p.Results.BandUpper)], ...
                    line_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                    'DisplayName', [p.Results.DisplayName, ' Confidence Band']);
            end
            
            % Hold on for additional plots
            hold on;
            grid on;
        end
        
        function cap_errors = capErrorValues(~, errors, max_value)
            % CAPERRORVALUES Cap error values to a maximum value
            %   cap_errors = CAPERRORVALUES(errors, max_value) caps error values in the 
            %   input to the specified maximum value
            %
            %   Input:
            %     errors: Input error values
            %     max_value: Maximum allowed value
            %   Output:
            %     cap_errors: Capped error values
            
            if iscell(errors)
                cap_errors = errors;
                [rows, cols] = size(errors);
                
                for i = 1:rows
                    for j = 1:cols
                        invalid_idx = isnan(errors{i,j}) | isinf(errors{i,j}) | (errors{i,j} > max_value);
                        cap_errors{i,j}(invalid_idx) = max_value;
                    end
                end
            else
                cap_errors = errors;
                invalid_idx = isnan(errors) | isinf(errors) | (errors > max_value);
                cap_errors(invalid_idx) = max_value;
            end
        end
    end
end
