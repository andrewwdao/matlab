classdef Metric < handle
    % Metric - A class for calculating various error metrics in statistical analysis
    % 
    % This class provides methods for calculating common error metrics used in 
    % performance evaluation of estimation algorithms, such as MSE, RMSE, median 
    % error, and percentile-based metrics.
    
    properties
        percentile_values = 50;  % Default percentiles to calculate - can also be a vector of 3 values [25, 50, 75]
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
        
        function plots(~, x_data, y_data, plot_type, varargin)
            % PLOTERRORMETRICS Create plot of error metrics with multiple lines and annotations
            %   PLOTERRORMETRICS(x_data, y_data, plot_type) creates a plot with y_data vs x_data
            %
            %   For single line:
            %     x_data: Vector of x-axis values
            %     y_data: Vector of y-axis values (error metrics) or matrix where each column is a line
            %
            %   plot_type: Type of plot ('linear', 'semilogy', 'loglog', etc.)
            %
            %   Optional Name-Value Pairs:
            %     'LineStyles' - Cell array of line styles (e.g., {'-', '--', ':', '-.'})
            %     'Markers' - Cell array of markers (e.g., {'o', 's', 'd', '^'})
            %     'LineWidth' - Line width (default: 2)
            %     'MarkerSize' - Marker size (default: 6)
            %     'DisplayNames' - Cell array of display names for legend
            %     'ShowBands' - Logical array indicating which lines should have error bands
            %     'BandLower' - Matrix where each column is lower band values for a line
            %     'BandUpper' - Matrix where each column is upper band values for a line
            %     'Title' - Title of the plot
            %     'XLabel' - X-axis label
            %     'YLabel' - Y-axis label
            %     'LegendLocation' - Location for the legend (default: 'northeast')
            %     'ShowAnnotation' - Whether to show annotation box (default: false)
            %     'AnnotationPosition' - Position [x y w h] (default: [0.15, 0.1, 0.3, 0.2])
            %     'AnnotationStrings' - Cell array of strings for annotation
            
            % Parse inputs
            p = inputParser;
            addParameter(p, 'LineStyles', {'-', '--', ':', '-.', '-', '--'});
            addParameter(p, 'Markers', {'o', 's', 'd', '^', 'v', 'p'});
            addParameter(p, 'LineWidth', 2);
            addParameter(p, 'MarkerSize', 6);
            addParameter(p, 'DisplayNames', {});
            addParameter(p, 'ShowBands', false);
            addParameter(p, 'BandLower', []);
            addParameter(p, 'BandUpper', []);
            addParameter(p, 'Title', '');
            addParameter(p, 'XLabel', 'Signal to Noise Ratio (SNR) [dB]');
            addParameter(p, 'YLabel', 'Error [m]');
            addParameter(p, 'LegendLocation', 'northeast');
            addParameter(p, 'ShowAnnotation', false);
            addParameter(p, 'AnnotationPosition', [0.14, 0.1, 0.3, 0.2]);
            addParameter(p, 'AnnotationStrings', {});
            parse(p, varargin{:});
            
            figure('Name', p.Results.Title, 'NumberTitle', 'off');

            % Get parameters
            lineStyles = p.Results.LineStyles;
            markers = p.Results.Markers;
            lineWidth = p.Results.LineWidth;
            markerSize = p.Results.MarkerSize;
            displayNames = p.Results.DisplayNames;
            showBands = p.Results.ShowBands;
            bandLower = p.Results.BandLower;
            bandUpper = p.Results.BandUpper;
            
            % Determine number of lines
            if isvector(y_data)
                numLines = 1;
                y_data = y_data(:); % Ensure column vector
            else
                [~, numLines] = size(y_data);
            end
            
            % Make sure x_data is compatible with y_data
            if isvector(x_data) && length(x_data) == size(y_data, 1)
                x_data = repmat(x_data(:), 1, numLines);
            end
            
            % Make sure displayNames is a cell array with enough elements
            if isempty(displayNames)
                displayNames = cell(1, numLines);
                for i = 1:numLines
                    displayNames{i} = ['Line ' num2str(i)];
                end
            elseif ~iscell(displayNames)
                displayNames = {displayNames};
            end
            
            % Extend single values to arrays if needed
            if isscalar(showBands)
                showBands = repmat(showBands, 1, numLines);
            end
            
            % Plot each line
            h_lines = zeros(numLines, 1);
            for i = 1:numLines
                % Select line style and marker
                lineStyle = lineStyles{mod(i-1, length(lineStyles))+1};
                marker = markers{mod(i-1, length(markers))+1};
                
                % Create plot based on type
                switch lower(plot_type)
                    case 'linear'
                        h_lines(i) = plot(x_data(:,i), y_data(:,i), ...
                            'LineStyle', lineStyle, ...
                            'Marker', marker, ...
                            'LineWidth', lineWidth, ...
                            'MarkerSize', markerSize, ...
                            'DisplayName', displayNames{i});
                    case 'semilogy'
                        h_lines(i) = semilogy(x_data(:,i), y_data(:,i), ...
                            'LineStyle', lineStyle, ...
                            'Marker', marker, ...
                            'LineWidth', lineWidth, ...
                            'MarkerSize', markerSize, ...
                            'DisplayName', displayNames{i});
                    case 'loglog'
                        h_lines(i) = loglog(x_data(:,i), y_data(:,i), ...
                            'LineStyle', lineStyle, ...
                            'Marker', marker, ...
                            'LineWidth', lineWidth, ...
                            'MarkerSize', markerSize, ...
                            'DisplayName', displayNames{i});
                    otherwise
                        error('Unknown plot type: %s', plot_type);
                end
                
                hold on;
                
                % Add confidence band if requested
                if showBands(i) && ~isempty(bandLower) && ~isempty(bandUpper)
                    % Get line color
                    lineColor = get(h_lines(i), 'Color');
                    
                    % Create shaded area
                    fill([x_data(:,i); flipud(x_data(:,i))], ...
                        [bandLower(:,i); flipud(bandUpper(:,i))], ...
                        lineColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                        'DisplayName', [displayNames{i}, ' Confidence Band']);
                end
            end
            
            % Add title, labels, and legend
            if ~isempty(p.Results.Title)
                title(p.Results.Title);
            end
            xlabel(p.Results.XLabel);
            ylabel(p.Results.YLabel);
            legend('Location', p.Results.LegendLocation);
            grid on;
            
            % Add annotation if requested
            if p.Results.ShowAnnotation && ~isempty(p.Results.AnnotationStrings)
                % Get annotation position and adjust based on number of lines
                annotPos = p.Results.AnnotationPosition;
                
                % Calculate required height based on number of annotation lines
                numLines = numel(p.Results.AnnotationStrings);
                baseHeight = 0.0325;  % Base height for box without text
                lineHeight = 0.0425; % Height per line of text
                requiredHeight = baseHeight + lineHeight * numLines;
                
                % Adjust y-position to ensure annotation stays in figure
                if annotPos(2) + requiredHeight > 1
                    % If would exceed top of figure, move down
                    annotPos(2) = max(0.01, 1 - requiredHeight - 0.01);
                end
                
                % Ensure bottom is not cut off
                if annotPos(2) < 0.05
                    annotPos(2) = 0.05;
                end
                
                % Update height based on content
                annotPos(4) = requiredHeight;
                
                % Create annotation with adjusted position
                annotation('textbox', annotPos, ...
                    'String', p.Results.AnnotationStrings, ...
                    'FitBoxToText', 'on', ...
                    'BackgroundColor', 'white', ...
                    'EdgeColor', 'black');
            end
        end

        
        function [result, cnt_capped, cnt_total] = capArrayValues(~, x, max_val, INCLUDE_CAPPED)
            % capArrayValues Cap values that are invalid or exceed maximum
            % Caps values in the input array that are NaN, Inf, or greater than the
            % specified maximum value.
            %
            % Syntax:
            %   [result, cnt_capped, cnt_total] = capArrayValues(x, max_val)
            %
            % Inputs:
            %   x       - Input array of any dimension
            %   max_val - Maximum value threshold
            %
            % Outputs:
            %   result     - Array with capped values
            %   cnt_capped - Number of elements that were capped
            %   cnt_total  - Total number of elements in the input array
            %
            % Example INCLUDE_CAPPED=true:
            %   [capped_data, num_capped, total] = capArrayValues([1, 5, NaN, Inf], 3);
            %   Returns: [1, 3, 3, 3]
            %   num_capped = 2, total = 4
            % Example with INCLUDE_CAPPED=false:
            %   [capped_data, num_capped, total] = capArrayValues([1, 5, NaN, Inf], 3, false);
            %   Returns: [1]
            %   num_capped = 2, total = 4
            invalid_idx = isnan(x) | isinf(x) | (x > max_val);
            cnt_capped = sum(invalid_idx);
            cnt_total = numel(x);
            if INCLUDE_CAPPED
                % Keep all values but cap the invalid ones
                result = x;
                result(invalid_idx) = max_val;
            else
                % Exclude the invalid values completely
                result = x(~invalid_idx);
            end
        end
        
        function result = capErrorValues(obj, errors, max_value, INCLUDE_CAPPED)
            % CAPERRORVALUES Cap error values and track capping statistics per method
            %   result = CAPERRORVALUES(errors, max_value) caps error values and returns
            %   statistics organized by method
            %
            %   Input:
            %     errors: Input error values (cell array or array)
            %     max_value: Maximum allowed value
            %   Output:
            %     result: Structure containing:
            %       .values - Capped error values (same structure as input)
            %       .cnt_capped - Cell array of capped counts per method or scalar value
            %       .cnt_total - Cell array of total counts per method or scalar value
            %       .percentage - Cell array of percentages per method or scalar value
            
            % Default behavior: include capped values
            if nargin < 4
                INCLUDE_CAPPED = true;
            end
            
            if iscell(errors)
                [rows, cols] = size(errors);
                
                % Process each cell element
                [errors_capped, capped_counts, total_counts] = cellfun(@(x) obj.capArrayValues(x, max_value, INCLUDE_CAPPED), ...
                    errors, 'UniformOutput', false);

                % Convert cell outputs to numeric matrices
                cnt_capped_detail = reshape(cell2mat(capped_counts), rows, cols);
                cnt_total_detail = reshape(cell2mat(total_counts), rows, cols);
                
                % Vectorized method-specific totals
                capped_sums = sum(cnt_capped_detail, 1);  % Sum each column
                total_sums = sum(cnt_total_detail, 1);    % Sum each column
                percentages = (capped_sums ./ total_sums) * 100;  % Calculate all percentages at once

                % Convert to cell arrays for output and Assign results
                result.cnt_capped = num2cell(capped_sums);
                result.cnt_total = num2cell(total_sums);
                result.percentage = num2cell(percentages);
                result.values = errors_capped;
            else
                % Handle non-cell input (standard behavior)
                invalid_idx = isnan(errors) | isinf(errors) | (errors > max_value);

                if INCLUDE_CAPPED
                    % Keep all values but cap the invalid ones
                    errors_capped = errors;
                    errors_capped(invalid_idx) = max_value;
                else
                    % Exclude the invalid values completely
                    errors_capped = errors(~invalid_idx);
                end

                % Track statistics for array input
                result.values = errors_capped;
                result.cnt_capped = sum(invalid_idx(:));
                result.cnt_total = numel(errors);
                result.percentage = (result.cnt_capped / result.cnt_total) * 100;
            end
        end
    end
end
