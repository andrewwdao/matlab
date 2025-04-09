function timer_output = progressbar(option, varargin)
%PROGRESSBAR Prints a progress bar.
%Syntax:
%   PROGRESSBAR('reset', total_progress)
%   PROGRESSBAR('resettimer')
%   PROGRESSBAR('step'[, step])
%   PROGRESSBAR('end')
%   PROGRESSBAR('barwidth', bar_width)
%   PROGRESSBAR('displaymode', mode)
%   PROGRESSBAR('minimalupdateinterval', min_interval);

persistent total_progress;
persistent current_progress;
persistent last_text_width;
persistent start_time;
persistent total_time;
persistent bar_width;
persistent display_mode;
persistent min_interval;
persistent last_print_time;
persistent timers;  % Structure to hold algorithm timers
persistent current_algo_name;  % Name of the currently running algorithm timer

if isempty(total_progress)
    total_progress = 0;
    current_progress = 0;
    last_text_width = 0;
    bar_width = 20;
    display_mode = 1;
    min_interval = 0.3;
    timers = struct();  % Initialize empty structure for algorithm timers
    current_algo_name = '';
end

switch lower(option)
    case 'minimalupdateinterval'
        min_interval = max(double(varargin{1}), 0);
    case 'displaymode'
        if nargin == 2
            switch lower(varargin{1})
                case 'replace'
                    display_mode = 0;
                case 'append'
                    display_mode = 1;
                otherwise
                    error('Display mode should be either ''replace'' (default) or ''append''.');
            end
        else
            error('Display mode not specified.');
        end
    case 'barwidth'
        if nargin == 2
            if ischar(varargin{1}) && strcmpi(varargin{1}, 'default')
                bar_width = 20;
            elseif isreal(varargin{1}) && varargin{1} > 0 && floor(varargin{1}) == varargin{1}
                bar_width = varargin{1};
            else
                error('Bar width must be a positive integer or ''default''.');
            end
        else
            error('Bar width is not specified.');
        end
    case 'resettimer'
        start_time = datetime();
        total_time = [];
    case 'reset'
        start_time = datetime();
        total_time = [];
        current_progress = 0;
        last_text_width = 0;
        if nargin == 2
            total_progress = varargin{1};
            if total_progress <= 0 || floor(total_progress) ~= total_progress
                error('Positive integer expected for the total progress.');
            end
        end
    case 'step'
        if nargin >= 2
            step = varargin{1};
            if isscalar(step) && isreal(step) && floor(step) == step
                current_progress = current_progress + step;
            else
                error('Step size should be an integer.');
            end
        else
            current_progress = current_progress + 1;
        end
        if current_progress > total_progress
            current_progress = total_progress;
        end
        now_time = datetime();
        % check if this function is called too frequently
        if current_progress < total_progress && ...
                ~isempty(last_print_time) && ...
                seconds(now_time - last_print_time) < min_interval
            return;
        end
        last_print_time = now_time;
        % print progress
        percentage = current_progress / total_progress;
        if percentage == 1.0        
            % completed
            if isempty(total_time)
                total_time = now_time - start_time;
            end
            text2print = sprintf('%s (%d/%d) %.1f%%%%\n[Runtime: %s]\n', ...
                get_progress_bar(percentage, bar_width), ...
                current_progress, total_progress, ...
                percentage * 100, seconds2str(seconds(total_time)));
        else
            % incomplete, generate ETA string
            if last_text_width == 0
                eta_str = '-';
            else
                seconds_elapsed = seconds(now_time - start_time);
                seconds_remaining = seconds_elapsed * ((1.0 - percentage) / percentage);
                eta_str = seconds2str(seconds_remaining);
            end
            text2print = sprintf('%s (%d/%d) %.1f%%%% ETA: %s\n', ...
                get_progress_bar(percentage, bar_width), ...
                current_progress, total_progress, ...
                percentage * 100, ...
                eta_str);
        end
        
        if last_text_width > 0 && display_mode == 0
            % not first print
            fprintf(repmat('\b', 1, last_text_width - 1));
        end
        fprintf(text2print);
        last_text_width = length(text2print);
    
    case 'starttimer'
        % Start timing an algorithm
        if nargin < 2
            error('Algorithm name must be provided');
        end
        
        algo_name = varargin{1};
        current_algo_name = algo_name;
        
        % Initialize timer for this algorithm if it doesn't exist
        if ~isfield(timers, algo_name)
            timers.(algo_name) = struct('times', [], 'current_start', []);
        end
        
        % Record start time
        timers.(algo_name).current_start = tic;
        
    case 'stoptimer'
        % Stop timing an algorithm and record the elapsed time
        if nargin < 2
            if isempty(current_algo_name)
                error('No algorithm timer is currently running');
            end
            algo_name = current_algo_name;
        else
            algo_name = varargin{1};
        end
        
        if ~isfield(timers, algo_name) || isempty(timers.(algo_name).current_start)
            error('Timer for algorithm "%s" was not started', algo_name);
        end
        
        % Calculate elapsed time
        elapsed = toc(timers.(algo_name).current_start);
        
        % Add to times array
        timers.(algo_name).times(end+1) = elapsed;
        timers.(algo_name).current_start = [];
        current_algo_name = '';
    
    case 'reporttimers'
        % Report average execution times for all algorithms
        fprintf('\n===== Algorithm Execution Times =====\n');
        
        algo_names = fieldnames(timers);
        if isempty(algo_names)
            fprintf('No algorithm timers have been recorded.\n');
        else
            for i = 1:length(algo_names)
                name = algo_names{i};
                times = timers.(name).times;
                
                if ~isempty(times)
                    avg_time = mean(times);
                    std_time = std(times);
                    fprintf('%s: %.4f s (±%.4f s) over %d runs\n', ...
                        name, avg_time, std_time, length(times));
                end
            end
        end
        fprintf('====================================\n');
    % Assign output
    timer_output = timers;
    
    case 'end'
        % force finish
        if current_progress < total_progress
            current_progress = total_progress;
            if isempty(start_time)
                error('End is called before starting.')
            end
            total_time = datetime() - start_time;
            if last_text_width > 0 && display_mode == 0
                % not first print
                fprintf(repmat('\b', 1, last_text_width - 1));
            end
            fprintf('%s (%d/%d) %.1f%%%%\n[Runtime: %s]\n', ...
                get_progress_bar(1.0, bar_width), ...
                current_progress, total_progress, ...
                100, seconds2str(seconds(total_time)));
        else
            % just print a new line
            fprintf('\n');
        end
    otherwise
        error('Unknown option.');
end

end

function bar = get_progress_bar(percent, width)
n_complete = floor(width * percent);
n_incomplete = width - n_complete;
if n_incomplete == 0
    bar = ['[' repmat('=', 1, n_complete) ']'];
else
    bar = ['[' repmat('=', 1, n_complete) '>' repmat(' ', 1, n_incomplete-1) ']'];
end
end

function s = seconds2str(t)
    %SECONDS2STR Convert seconds to nicer looking string.
    t = double(t);
    if t < 1e-3
        s = sprintf('%3.2f ns', t*1e6);
    elseif t < 1
        s = sprintf('%3.2f ms', t*1e3);
    elseif t < 60
        s = sprintf('%2.2f s', t);
    else
        if t > 31536000
            year = floor(t/31536000);
            t = mod(t, 31536000);
        else
            year = 0;
        end
        if t > 86400
            day = floor(t/86400);
            t = mod(t, 86400);
        else
            day = 0;
        end
        if t > 3600
            hour = floor(t/3600);
            t = mod(t, 3600);
        else
            hour = 0;
        end
        minute = floor(t/60);
        second = floor(mod(t, 60));
        s = sprintf('%d m %d s', minute, second);
        if hour > 0
            s = [sprintf('%d h ', hour) s];
        end
        if day > 0
            s = [sprintf('%d d ', day) s];
        end
        if year > 0
            s = [sprintf('%d y ', day) s];
        end
    end
    
end
