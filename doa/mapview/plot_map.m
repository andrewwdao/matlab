% Define the first point and angle
x0_1 = 1; % x-coordinate of the first point
y0_1 = 2; % y-coordinate of the first point
angle1 = 45; % angle in degrees with respect to the x-axis

% Define the second point and angle
x0_2 = 4; % x-coordinate of the second point
y0_2 = 1; % y-coordinate of the second point
angle2 = 135; % angle in degrees with respect to the x-axis

% Convert angles to radians
theta1 = deg2rad(angle1);
theta2 = deg2rad(angle2);

% Define the line equations: y = mx + c
m1 = tan(theta1);
c1 = y0_1 - m1 * x0_1;

m2 = tan(theta2);
c2 = y0_2 - m2 * x0_2;

% Calculate intersection point
xi = (c2 - c1) / (m1 - m2);
yi = m1 * xi + c1;

% Plot the lines
figure;
hold on;
fplot(@(x) m1*x + c1, [x0_1 10], 'r'); % First line
fplot(@(x) m2*x + c2, [-4 x0_2], 'b'); % Second line
plot(xi, yi, 'ko', 'MarkerSize', 10, 'LineWidth', 2); % Intersection point
xlabel('X');
ylabel('Y');
legend('Line 1', 'Line 2', 'Intersection Point');
title('Intersection of Two Infinite Lines');
grid on;
axis equal;
hold off;

% Display the intersection point
disp(['Intersection Point: (', num2str(xi), ', ', num2str(yi), ')']);
