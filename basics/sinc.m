% MATLAB Code to Plot the Sinc Function with Only x=0 and y=0 Axes

% Define the range of x values
x = linspace(-10, 10, 1000); % 1000 points between -10 and 10

% Compute the sinc function
y = sinc(x); % MATLAB has a built-in sinc function

% Plot the sinc function
figure; % Create a new figure
plot(x, y, 'k-', 'LineWidth', 1.5); % Plot with a blue line
hold on; % Hold the plot for additional elements

% Add x=0 and y=0 axes
plot([-10, 10], [0, 0], 'k-', 'LineWidth', 1.2); % Draw the x-axis (y=0)
plot([0, 0], [-0.3, 1.1], 'k-', 'LineWidth', 1.2); % Draw the y-axis (x=0)

% Customize the plot
grid off; % Turn off the grid
axis off; % Turn off the default axis box
xlim([-10, 10]); % Set x-axis limits
ylim([-0.3, 1.1]); % Set y-axis limits

% Add labels for x=0 and y=0 axes
text(10.2, 0, '\tau', 'FontSize', 12, 'FontWeight', 'bold'); % Label for x-axis
text(0, 1.15, 'R_n(\tau)', 'FontSize', 12, 'FontWeight', 'bold'); % Label for y-axis

% title('Sinc Function with Only x=0 and y=0 Axes'); % Add a title