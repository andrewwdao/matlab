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
            %  Calculate Mean Square Error
            %  from the input error values
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
                        valid_errors = errors{i, j}(~isnan(errors{i, j}) & ~isinf(errors{i, j}));
                        if isempty(valid_errors)
                            mse(i, j) = NaN; % If no valid values remain, result is NaN
                        else
                            % Handle potential overflow during squaring
                            max_safe_value = sqrt(realmax);
                            if any(abs(valid_errors) > max_safe_value)
                                % Option 1: Cap extremely large values
                                too_large = abs(valid_errors) > max_safe_value;
                                valid_errors(too_large) = sign(valid_errors(too_large)) * max_safe_value;
                                warning('Overflow detected in squaring errors. Values capped to avoid overflow.');
                            end
                            mse(i, j) = mean(valid_errors.^2);
                        end
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
                
                % Initialize structure for 3 percentiles, assume [lower, median, upper]
                percentiles = struct('lower', zeros(rows, cols), ...
                                    'val', zeros(rows, cols), ...
                                    'upper', zeros(rows, cols));
                % Sort percentile values to ensure lower/middle/upper order
                percentile_values = sort(percentile_values);
                
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
            addParameter(p, 'Title', 'Error Metric by estimation method');
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
                numLines = size(y_data, 2);
            end
            fprintf('Number of lines to plot: %d\n', numLines);
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
                        fprintf('Plotting semilogy for line %d\n', i);
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
                        'DisplayName', [displayNames{i}, ' Variance']);
                end
            end
            
            % Add title, labels, and legend
            title(p.Results.Title);
            xlabel(p.Results.XLabel); xlim([min(min(x_data)), max(max(x_data))]);
            ylabel(p.Results.YLabel);
            legend('Location', p.Results.LegendLocation);
            grid on; xticks(unique(x_data(:))); % Add x-axis grid lines at each x_data point
            
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

        function plot_epdf(obj, all_errors, legend_names, snr_indices, figure_title, varargin)
            % PLOT_EPDF Plot empirical PDF of position errors for different methods
            %
            % Visualization of error distributions with additional statistical metrics
            % for comparing estimation methods.
            %
            % Parameters:
            %   all_errors - Cell array of error values from Monte Carlo simulation
            %   legend_names - Cell array of method names
            %   snr_indices - SNR indices to plot (default = [1, ceil(end/2), end])
            %   figure_title - Title for the figure (optional)
            %   varargin - Additional options as name-value pairs:
            %       'BinWidth' - Width of histogram bins (default = auto)
            %       'CapInfo' - Optional struct with capping information
            %       'SNRValues' - SNR values in dB corresponding to indices
            %       'ShowKDE' - Whether to show kernel density estimation (default = true)
            %       'CompareMethod' - Index of reference method for comparison (default = [])
            %       'CompareInSubplot' - Add comparison as subplot (true) or separate figure (false) (default = false)
            %       'ColorMap' - Colormap for plots (default = 'parula')
            
            % Parse inputs
            p = inputParser;
            addParameter(p, 'BinWidth', []);  % Auto-calculate bin width if not specified
            addParameter(p, 'CapInfo', []);
            addParameter(p, 'SNRValues', []);
            addParameter(p, 'ShowKDE', true); % Show KDE curve by default
            addParameter(p, 'CompareMethod', []);
            addParameter(p, 'CompareInSubplot', false); % Default to separate figure
            addParameter(p, 'ColorMap', 'parula');
            parse(p, varargin{:});
            
            % Extract parameters
            cap_info = p.Results.CapInfo;
            snr_values = p.Results.SNRValues;
            show_kde = p.Results.ShowKDE;
            compare_method = p.Results.CompareMethod;
            compare_in_subplot = p.Results.CompareInSubplot;
            cmap = colormap(p.Results.ColorMap);
            
            % Default SNR indices if not provided
            if nargin < 4 || isempty(snr_indices)
            snr_indices = [1, ceil(size(all_errors, 1)/2), size(all_errors, 1)];
            end
            
            % Default figure title
            if nargin < 5 || isempty(figure_title)
            figure_title = 'Empirical PDF of Position Errors';
            end
            
            % Get number of methods and SNR levels to plot
            num_methods = size(all_errors, 2);
            num_snr_levels = length(snr_indices);
            
            % Create a figure for each SNR level
            for snr_idx_pos = 1:num_snr_levels
            snr_idx = snr_indices(snr_idx_pos);
            
            % Create figure with title including SNR information
            if ~isempty(snr_values)
                fig_title = sprintf('%s (SNR: %.1f dB)', figure_title, snr_values(snr_idx));
            else
                fig_title = sprintf('%s (SNR Index: %d)', figure_title, snr_idx);
            end
            
            % Create figure with optimized size
            figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
            
            % Calculate optimal subplot layout
            % If we're including a comparison subplot, add one to the total count
            total_subplots = num_methods;
            if ~isempty(compare_method) && compare_in_subplot
                total_subplots = total_subplots + 1;
            end
            [rows, cols] = obj.subplot_layout(total_subplots);
            
            global_max_error = 0;
            % Create subplots for each method
            for method_idx = 1:num_methods
                % Extract errors for current method at specified SNR
                errors = all_errors{snr_idx, method_idx};
                method_max_error = min(max(errors), 1.5*mean(errors(:)));
                % Auto-calculate bin width if not specified
                if isempty(p.Results.BinWidth)
                    % Calculate method-specific bin width based on its error range
                    method_bin_width = (method_max_error/30); % 30 bins for each method
                    bin_width = method_bin_width;
                else
                    bin_width = p.Results.BinWidth;
                end
                % Calculate statistics
                mean_error = mean(errors);
                median_error = median(errors);
                std_error = std(errors);
                rmse = sqrt(mean(errors.^2));
                iqr_value = iqr(errors);
                
                % Create subplot
                subplot(rows, cols, method_idx);
                
                % Determine color for current method
                color_idx = mod(method_idx-1, size(cmap, 1)) + 1;
                method_color = cmap(color_idx, :);
                
                % Plot histogram
                histogram(errors, 'BinWidth', bin_width, ...
                'Normalization', 'pdf', ...
                'FaceColor', method_color, ...
                'FaceAlpha', 0.7, ...
                'EdgeColor', 'k');
                
                hold on;
                
                % Add KDE curve if requested
                % if show_kde
                %     % Add a small epsilon to zero values
                %     epsilon = 1e-10;
                %     positive_errors = errors;
                %     positive_errors(positive_errors == 0) = epsilon;
                    
                %     try
                %         if ~isempty(positive_errors)
                %             % Analyze distribution characteristics
                %             error_skewness = skewness(positive_errors);
                %             error_kurtosis = kurtosis(positive_errors);
                %             cv = std(positive_errors) / mean(positive_errors); % coefficient of variation
                            
                %             % Select appropriate KDE method based on distribution characteristics
                %             if cv > 1.2 || error_skewness > 2
                %                 % Heavy-tailed or highly skewed data - use log-transform
                %                 log_errors = log(positive_errors);
                %                 [f_log, xi_log] = ksdensity(log_errors, 'Kernel', 'epanechnikov');
                                
                %                 % Transform back to original scale
                %                 xi = exp(xi_log);
                %                 f = f_log ./ xi; % Adjustment for variable transformation
                %             elseif error_kurtosis > 3 && error_kurtosis <= 5
                %                 % Moderately heavy-tailed - use Epanechnikov kernel with larger bandwidth
                %                 bw = 1.2 * std(positive_errors) * length(positive_errors)^(-0.2);
                %                 [f, xi] = ksdensity(positive_errors, 'Kernel', 'epanechnikov', 'Bandwidth', bw, 'Support', 'positive');
                %             elseif length(positive_errors) < 100
                %                 % Small sample size - use robust bandwidth
                %                 bw = 1.06 * min(std(positive_errors), iqr(positive_errors)/1.34) * length(positive_errors)^(-0.2);
                %                 [f, xi] = ksdensity(positive_errors, 'Bandwidth', bw, 'Support', 'positive');
                %             else
                %                 % Default case - adaptive bandwidth calculation manually
                %                 bw_silverman = 0.9 * min(std(positive_errors), iqr(positive_errors)/1.34) * length(positive_errors)^(-0.2);
                %                 [f, xi] = ksdensity(positive_errors, 'Bandwidth', bw_silverman, 'Support', 'positive');
                %             end
                            
                %             % Scale the KDE curve to match histogram height
                %             hist_counts = histcounts(positive_errors, 'Normalization', 'pdf');
                %             if ~isempty(hist_counts)
                %                 histMax = max(hist_counts);
                %                 kdeMax = max(f);
                %                 if kdeMax > 0
                %                     scaling_factor = histMax / kdeMax;
                %                     f = f * scaling_factor;
                %                 end
                %             end
                            
                %             % Fit theoretical distributions and find best match
                %             [best_dist, best_params, gof] = obj.findBestDistribution(positive_errors);
                            
                %             % Plot KDE curve
                %             plot(xi, f, 'LineWidth', 2, 'Color', method_color*0.7); hold on;
                            
                %             % Add distribution name to title or as text annotation
                %             if exist('title_text', 'var')
                %                 % Format parameters differently based on distribution type
                %                 param_display = obj.formatDistParams(best_dist, best_params);
                %                 title_text = sprintf('%s\nBest fit: %s %s (GoF: %.3f)', title_text, best_dist, param_display, gof);
                %                 title(title_text);
                %             else
                %                 % Format parameters for annotation
                %                 param_display = obj.formatDistParams(best_dist, best_params);
                %                 text(0.5*max(xi), 0.7*max(f), sprintf('Best fit: %s %s\nGoF: %.3f', best_dist, param_display, gof), 'FontSize', 8);
                                
                %                 % Optionally plot the theoretical distribution
                %                 x_range = linspace(0, max(xi), 200);
                %                 if ~isempty(best_params) && ~strcmp(best_dist, 'unknown')
                %                     y_theo = obj.distPDF(x_range, best_dist, best_params);
                %                     % Scale theoretical PDF to match histogram height
                %                     y_theo = y_theo * (max(f)/max(y_theo));
                %                     plot(x_range, y_theo, '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
                %                 end
                %             end
                %             % Print distribution information to the terminal
                %             method_name = '';
                %             if exist('legend_names', 'var') && method_idx <= length(legend_names)
                %                 method_name = legend_names{method_idx};
                %             end

                %             % Format and print the distribution info to terminal
                %             param_display = obj.formatDistParams(best_dist, best_params);
                %             fprintf('Method: %s, SNR Index: %d\n', method_name, snr_idx);
                %             fprintf('  Best fit distribution: %s %s\n', best_dist, param_display);
                %             fprintf('  Goodness of fit: %.3f\n', gof);
                %             fprintf('  Distribution parameters: ');
                %             disp(best_params);
                %             fprintf('--------------------------\n');
                %         end
                %     catch e
                %         warning(e.identifier, 'KDE calculation failed: %s', e.message);
                %     end
                % end
                
                % Create detailed title with statistics
                if ~isempty(cap_info)
                cap_percentage = cap_info.percentage{method_idx};
                title_text = sprintf('%s\nMean: %.1f m, Std: %.1f m\nMedian: %.1f m, RMSE: %.1f m\nIQR: %.1f m, Capped: %.1f%%', ...
                    legend_names{method_idx}, mean_error, std_error, median_error, rmse, iqr_value, cap_percentage);
                else
                title_text = sprintf('%s\nMean: %.1f m, Std: %.1f m\nMedian: %.1f m, RMSE: %.1f m\nIQR: %.1f m', ...
                    legend_names{method_idx}, mean_error, std_error, median_error, rmse, iqr_value);
                end
                title(title_text);
                
                % Labels
                xlabel('Position Error (m)');
                ylabel('Probability Density');
                grid on;
                
                % Set consistent axis limits
                xlim([0, method_max_error]);
                global_max_error = max(global_max_error, method_max_error);
            end
            
            % Add comparison plot 
            if ~isempty(compare_method) && compare_method <= num_methods
                if compare_in_subplot
                % Add comparison as final subplot
                subplot(rows, cols, num_methods + 1);
                obj.plot_comparison_subplot(all_errors, snr_idx, compare_method, legend_names, global_max_error, snr_values);
                else
                % Create a separate figure for method comparison
                figure('Name', [fig_title, ' - Method Comparison'], 'NumberTitle', 'off', 'Position', [150, 150, 1000, 600]);
                obj.plot_comparison_subplot(all_errors, snr_idx, compare_method, legend_names, global_max_error, snr_values);
                end
            end
            
            % Add overall title to main figure
            sgtitle(fig_title, 'FontSize', 14, 'FontWeight', 'bold');
            end
        end

        % Function to find the best theoretical distribution
        function [best_dist, best_params, best_gof] = findBestDistribution(~, data)
            % List of distributions to try
            dists = {'gamma', 'rayleigh', 'weibull', 'lognormal', 'exponential'};
            
            best_gof = Inf;
            best_dist = '';
            best_params = [];
            
            for i = 1:length(dists)
                dist_name = dists{i};
                try
                    % Fit distribution
                    pd = fitdist(data, dist_name);
                    
                    % Calculate goodness of fit using Kolmogorov-Smirnov test
                    [~, p] = kstest(data, 'CDF', [data, cdf(pd, data)]);
                    
                    % Alternative: use negative log-likelihood as goodness of fit
                    % nll = negloglik(pd);
                    
                    % Calculate AIC (Akaike Information Criterion)
                    if p < 0.05
                        % Poor fit, higher AIC penalty
                        aic = pd.NLogL + 4*length(pd.Params);
                    else
                        % Good fit by KS test
                        aic = pd.NLogL + 2*length(pd.Params);
                    end
                    
                    % If this distribution has better fit (smaller AIC)
                    if aic < best_gof
                        best_gof = aic;
                        best_dist = dist_name;
                        best_params = pd.Params;
                    end
                catch
                    % Skip distributions that fail to fit
                    continue;
                end
            end
            
            % If no distribution fit worked, return default
            if isempty(best_dist)
                best_dist = 'unknown';
                best_params = [];
                best_gof = Inf;
            end
        end

        function param_display = formatDistParams(~, dist_name, params)
            % Format distribution parameters for display
            if isempty(params)
                param_display = '';
                return;
            end
            
            switch dist_name
                case 'gamma'
                    param_display = sprintf('(α=%.2f, β=%.2f)', params(1), params(2));
                case 'rayleigh'
                    param_display = sprintf('(σ=%.2f)', params);
                case 'weibull'
                    param_display = sprintf('(α=%.2f, β=%.2f)', params(1), params(2));
                case 'lognormal'
                    param_display = sprintf('(μ=%.2f, σ=%.2f)', params(1), params(2));
                case 'exponential'
                    param_display = sprintf('(λ=%.2f)', params);
                otherwise
                    param_display = '';
            end
        end
        
        function y = distPDF(~, x, dist_name, params)
            % Generate PDF values for a given distribution
            switch dist_name
                case 'gamma'
                    y = gampdf(x, params(1), params(2));
                case 'rayleigh'
                    y = raylpdf(x, params);
                case 'weibull'
                    y = wblpdf(x, params(1), params(2));
                case 'lognormal'
                    y = lognpdf(x, params(1), params(2));
                case 'exponential'
                    y = exppdf(x, 1/params); % Convert rate to mean
                otherwise
                    y = zeros(size(x));
            end
        end
        
        function plot_comparison_subplot(~, all_errors, snr_idx, compare_method, legend_names, method_max_error, snr_values)
            % Helper function to create comparison plot (either as subplot or separate figure)
            
            % Get reference errors
            ref_errors = all_errors{snr_idx, compare_method};
            
            % Calculate KDE for reference method
            % Add a small epsilon to zero values
            epsilon = 1e-10; % Small enough to be negligible for your position error scale
            adjusted_ref_errors = ref_errors;
            adjusted_ref_errors(adjusted_ref_errors == 0) = epsilon;
            
            % Now all values are positive
            [f_ref, xi_ref] = ksdensity(adjusted_ref_errors, 'Support', 'positive');
            
            % Plot reference KDE
            plot(xi_ref, f_ref, 'LineWidth', 2, 'Color', 'k', 'DisplayName', legend_names{compare_method}); 
            hold on;
            
            % Plot reference mean line
            ref_mean = mean(ref_errors);
            xline(ref_mean, 'k--', '', 'LineWidth', 1.5);
            
            % Initialize labels
            plot_labels = cell(1, size(all_errors, 2)+1);
            plot_labels{1} = legend_names{compare_method};
            plot_labels{2} = sprintf('Mean %s', legend_names{compare_method});
            
            % Track position in labels array
            label_idx = 3;
            
            % Plot KDEs for all other methods
            for method_idx = 1:size(all_errors, 2)
            % Skip reference method
            if method_idx == compare_method
                continue;
            end
            
            % Get errors and calculate KDE
            % Add a small epsilon to zero values
            epsilon = 1e-10; % Small enough to be negligible for your position error scale
            curr_errors = all_errors{snr_idx, method_idx};
            curr_errors(curr_errors == 0) = epsilon;
            
            % Now all values are positive
            [f_curr, xi_curr] = ksdensity(curr_errors, 'Support', 'positive');
            % Plot KDE
            plot(xi_curr, f_curr, 'LineWidth', 2);
            
            % Update labels
            plot_labels{label_idx} = legend_names{method_idx};
            label_idx = label_idx + 1;
            end
            
            % Clean up plot
            grid on;
            xlabel('Position Error (m)');
            ylabel('Probability Density');
            
            % Remove empty elements from labels
            plot_labels = plot_labels(~cellfun('isempty', plot_labels));
            
            % Add legend
            legend(plot_labels, 'Location', 'best');
            
            % Add title
            if ~isempty(snr_values)
            title(sprintf('Error Distribution Comparison at SNR = %.1f dB', snr_values(snr_idx)));
            else
            title(sprintf('Error Distribution Comparison at SNR Index = %d', snr_idx));
            end
            
            % Set x-axis limits
            xlim([0, method_max_error]);
        end
        
        function [rows, cols] = subplot_layout(~, n)
            % Determine optimal subplot layout based on number of plots
            
            % Calculate possible arrangements
            sqrtn = sqrt(n);
            cols = ceil(sqrtn);
            rows = ceil(n/cols);
            
            % Ensure we have enough spaces
            if rows*cols < n
            cols = cols + 1;
            end
        end
        
        function tight_layout(~)
            % Apply tight layout to current figure
            set(gcf, 'Units', 'normalized');
            set(findall(gcf, 'Type', 'axes'), 'LooseInset', get(gca, 'TightInset'));
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
                % If all values were excluded, set result to 0
                % if isempty(result)
                %     result = max_val;  % Set to max value to indicate all values were capped
                % end
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
