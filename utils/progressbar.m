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
        timers = struct();  % Initialize empty structure for algorithm timers
        current_algo_name = '';
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
        % Start timing an algorithm with flexible structure
        % Format: 'parent' or 'parent/child' or 'parent/child/grandchild'
        if nargin < 2
            error('Algorithm name must be provided');
        end
        
        algo_name = varargin{1};
        current_algo_name = algo_name;
        
        % Parse structure using '/' as separator
        parts = strsplit(algo_name, '/');
        
        % Initialize parent if it doesn't exist
        parent = parts{1};
        if ~isfield(timers, parent)
            timers.(parent) = struct();
        end
        
        % Handle single-level, 2-level, or 3-level hierarchy
        if isscalar(parts)
            % Single-level: just parent
            if ~isfield(timers.(parent), 'times')
                timers.(parent).times = [];
                timers.(parent).current_start = [];
            end
            timers.(parent).current_start = tic;
        elseif length(parts) == 2
            % Two-level: parent/child
            child = parts{2};
            if ~isfield(timers.(parent), child)
                timers.(parent).(child) = struct('times', [], 'current_start', []);
            end
            timers.(parent).(child).current_start = tic;
        elseif length(parts) == 3
            % Three-level: parent/child/grandchild
            child = parts{2};
            grandchild = parts{3};
            if ~isfield(timers.(parent), child)
                timers.(parent).(child) = struct();
            end
            if ~isfield(timers.(parent).(child), grandchild)
                timers.(parent).(child).(grandchild) = struct('times', [], 'current_start', []);
            end
            timers.(parent).(child).(grandchild).current_start = tic;
        else
            error('Timer hierarchy must be either "parent", "parent/child", or "parent/child/grandchild"');
        end
    
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
        
        % Parse parent/child structure
        parts = strsplit(algo_name, '/');
        parent = parts{1};
        
        % Check if parent exists
        if ~isfield(timers, parent)
            error('Timer for algorithm "%s" was not started', algo_name);
        end
        
        % Handle single-level, 2-level, or 3-level hierarchy
        if isscalar(parts)
            % Single-level: just parent
            if ~isfield(timers.(parent), 'current_start') || isempty(timers.(parent).current_start)
                error('Timer for algorithm "%s" was not started', algo_name);
            end
            elapsed = toc(timers.(parent).current_start);
            timers.(parent).times(end+1) = elapsed;
            timers.(parent).current_start = [];
        elseif length(parts) == 2
            % Two-level: parent/child
            child = parts{2};
            if ~isfield(timers.(parent), child) || isempty(timers.(parent).(child).current_start)
                error('Timer for algorithm "%s" was not started', algo_name);
            end
            elapsed = toc(timers.(parent).(child).current_start);
            timers.(parent).(child).times(end+1) = elapsed;
            timers.(parent).(child).current_start = [];
        elseif length(parts) == 3
            % Three-level: parent/child/grandchild
            child = parts{2};
            grandchild = parts{3};
            if ~isfield(timers.(parent), child) || ~isfield(timers.(parent).(child), grandchild) || ...
            isempty(timers.(parent).(child).(grandchild).current_start)
                error('Timer for algorithm "%s" was not started', algo_name);
            end
            elapsed = toc(timers.(parent).(child).(grandchild).current_start);
            timers.(parent).(child).(grandchild).times(end+1) = elapsed;
            timers.(parent).(child).(grandchild).current_start = [];
        end
        
        current_algo_name = '';

    case 'reporttimers'
        % Report average execution times for all algorithms
        fprintf('\n===== Algorithm Execution Times =====\n');
        
        % Initialize flattened output structure array
        flat_timer_stats = struct('FullName', {}, 'Parent', {}, 'Name', {}, 'Mean', {}, 'Std', {}, 'Count', {}, 'Times', {});
        
        % Call recursive function to populate flat_timer_stats
        if ~isempty(timers) && isstruct(timers)
            timers_struct = timers; % Use a local copy
            flat_timer_stats = process_timer_node(timers_struct, '', flat_timer_stats);
        end

        % --- Print the report (using the flat structure for easier iteration) ---
        if isempty(flat_timer_stats)
            fprintf('No algorithm timers have been recorded.\n');
        else
            % Sort by FullName for consistent printing order (optional)
            try % Sorting might fail if FullName is not always a string
                [~, sort_idx] = sort({flat_timer_stats.FullName});
                flat_timer_stats = flat_timer_stats(sort_idx);
            catch
                % Proceed without sorting if error occurs
            end

            for i = 1:length(flat_timer_stats)
                stat = flat_timer_stats(i);
                % Ensure fields exist before accessing
                if isfield(stat, 'FullName') && isfield(stat, 'Mean') && isfield(stat, 'Std') && isfield(stat, 'Count')
                    display_name = strrep(stat.FullName, '_', ' '); % Use FullName for display
                    fprintf('%s: %.4f s (±%.4f s) over %d runs\n', ...
                        display_name, stat.Mean, stat.Std, stat.Count);
                end
            end
        end
        
        fprintf('====================================\n');
        
        % Assign the *flattened* structure array as the output
        timer_output = flat_timer_stats; % Return the flat structure
        
    
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

function flat_stats = process_timer_node(node, current_path, flat_stats)
    if ~isstruct(node)
        return; % Stop if node is not a structure
    end
    fields = fieldnames(node);
    for i = 1:length(fields)
        field_name = fields{i};
        % Skip internal fields if they exist at this level by mistake
        if strcmp(field_name, 'times') || strcmp(field_name, 'current_start')
            continue;
        end
        
        sub_node = node.(field_name);
        
        % Construct the full path for this node
        if isempty(current_path)
            full_name = field_name;
        else
            full_name = [current_path '/' field_name];
        end

        % Check if this is a leaf node (contains 'times' field)
        if isstruct(sub_node) && isfield(sub_node, 'times') && ~isempty(sub_node.times)
            times = sub_node.times;
            % Ensure times is numeric before calculations
            if isnumeric(times)
                avg_time = mean(times);
                std_time = std(times);
                count = length(times);
                
                % Extract parent and name
                path_parts = strsplit(full_name, '/');
                if length(path_parts) > 1
                    parent_name = strjoin(path_parts(1:end-1), '/');
                    leaf_name = path_parts{end};
                else
                    parent_name = ''; % Top-level timer
                    leaf_name = full_name;
                end

                % Add entry to the flat structure array
                new_entry = struct(...
                    'FullName', full_name, ...
                    'Parent', parent_name, ...
                    'Name', leaf_name, ...
                    'Mean', avg_time, ...
                    'Std', std_time, ...
                    'Count', count, ...
                    'Times', times);
                flat_stats(end+1) = new_entry;
            end
            
        elseif isstruct(sub_node) % It's an intermediate node, recurse
            flat_stats = process_timer_node(sub_node, full_name, flat_stats);
        end
    end
end
