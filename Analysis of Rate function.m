% Define the parameters
k1 = 1.0; k2 = 0.5; %assuming these parameters as not given in question
% Define the domain
Ca_strt = 0.0; Ca_end = 5.0; Ca_step = 0.01;
syms Ca
rate = @(Ca) (k1 * Ca)./ (1 + k2 * (Ca.^2));

% (a) Derivative at a point
Ca_pt = 2.0;
slopex =  ((rate(Ca_pt+0.01)-rate(Ca_pt))/0.01);
fprintf('(a) Slope of expression at Ca = 2.0 is %f\n', slopex);

%(b) Calculating the area under curve using trapz in-built function from
%0.0 to 5.0
n_points = 501; % number of points in the vector
Ca1 = linspace(Ca_strt, Ca_end, n_points); % create the vector
auc = trapz(Ca1, rate(Ca1));
fprintf('(b) Area under curve from 0.0 to 5.0 is %f\n', auc);

% (c) Finding a minimum or maximum functional value over the given domain
rate_min = min(rate(Ca1));
rate_max = max(rate(Ca1));
fprintf('(c) Minimum Value of expression is %f\n',rate_min);
fprintf('(c) Maximum Value of expression is %f\n',rate_max);

% (d) Drawing a straight line between two points on the function
Ca_start_line = 0.5; Ca_end_line = 1.5; % Define the points
rate_start_line = rate(Ca_start_line); rate_end_line = rate(Ca_end_line);
m = (rate_end_line - rate_start_line) / (Ca_end_line - Ca_start_line); % Calculate the slope
c = rate_start_line - m * Ca_start_line; % Calculate the intercept
line_x = [Ca_start_line, Ca_end_line];
line_y = m * line_x + c;
eqn = ['y = ', num2str(m), 'x + ', num2str(c)];
fprintf('(d) Equation of line of two given points\nx1 = %f, y1 = %f, x2 = %f, y2 = %f is\n', Ca_start_line, rate_start_line, Ca_end_line,rate_end_line); 
disp(eqn) %print the equation

% (e)   Searching for a point where the tangent is the same as the given slope
slope = m; % Define the slope
Ca_tangent = Ca_strt;
max_iter = 1000; % Define the maximum number of iterations
iter = 0;
epsilon = 0.001; % Define the convergence criterion
while abs(((rate(Ca_tangent+0.01)-rate(Ca_tangent))/0.01) - slope) > epsilon && iter < max_iter
    Ca_tangent = Ca_tangent + epsilon; %Update the Ca value 
    iter = iter + 1;
end
fprintf("(e) Point on the curve where tangent is same as given slope is\nx = %f, y = %f\n", Ca_tangent,rate(Ca_tangent));
tangent_slope = ((rate(Ca_tangent+0.01)-rate(Ca_tangent))/0.01);

% Find the y-intercept of the line with slope m passing through the point Ca_tangent, slope_function(Ca_tangent, k1, K2)
y1_intercept = rate(Ca_tangent) - m*Ca_tangent;

% Define the x values for the line
line_x1 = linspace(Ca_strt, Ca_end, 100);

% Define the y values for the line using y = mx + b
line_y1 = m*line_x1 + y1_intercept;

% Plot the rate expression curve and the line
figure();
plot(Ca1, rate(Ca1), 'b','LineWidth',1);
hold on;
plot(line_x, line_y, 'r');
plot(Ca_tangent, rate(Ca_tangent) , 'ro'); % Plot the point where the tangent is the same as the given slope
% Plot the rate expression curve and the line
plot(line_x1, line_y1, 'r');
xlabel('Ca');
ylabel('-rA');
legend('Rate expression', 'Line', 'Point where tangent = Given slope');
title('Rate expression curve and Line');