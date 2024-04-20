% Synthetic Data Generation
% Let's consider a simple quadratic equation: y = ax^2 + bx + c
a = 2;
b = -3;
c = 1;

% Generate synthetic data with noise
x = linspace(-10, 10, 100);
y_true = a * x.^2 + b * x + c;
noise = 5 * randn(size(x));
y_observed = y_true + noise;

% Plot the synthetic data
figure;
plot(x, y_observed, 'bo', 'MarkerFaceColor', 'b');
hold on;
plot(x, y_true, 'r', 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('Synthetic Data');
legend('Observed Data', 'True Model');

% Define Objective Function (Cost Function)
objective_function = @(params) sum((params(1) * x.^2 + params(2) * x + params(3) - y_observed).^2);

% Global Optimization (Simulated Annealing)
options = optimoptions(@simulannealbnd, 'Display', 'iter');
lb = [-10, -10, -10];  % Lower bounds for parameters
ub = [10, 10, 10];     % Upper bounds for parameters
initial_guess = [1, 1, 1];  % Initial guess for parameters

[params_optimized, fval] = simulannealbnd(objective_function, initial_guess, lb, ub, options);

% Plot the optimized model
y_optimized = params_optimized(1) * x.^2 + params_optimized(2) * x + params_optimized(3);
plot(x, y_optimized, 'g', 'LineWidth', 2);
legend('Observed Data', 'True Model', 'Optimized Model');
hold off;

disp('Optimized Parameters:');
disp(params_optimized);
disp(['Optimized Cost (RSS): ', num2str(fval)]);