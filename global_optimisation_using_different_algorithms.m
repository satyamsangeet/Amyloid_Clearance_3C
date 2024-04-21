%Defining example ODE, y = ax^2 + bx + c
a = 2;
b = -1;
c = 1;

%generating synthetic dataset
x = linspace(-10, 10, 100);
y_true = a * x.^2 + b * x + c;
noise = 5 * randn(size(x));
y_observed = y_true + noise;

%plotting synthetic data
figure;
plot(x, y_observed, 'bo', 'MarkerFaceColor', 'b');
hold on;
plot(x, y_true, 'r', 'LineWidth', 2);
xlabel('x');
ylabel('y');
legend('Observed Data', 'True Model');

%Defining objective function
objective_function = @(params) sum((params(1) * x.^2 + params(2) * x + params(3) - y_observed).^2);

%Global optimisation using Simulated Annealing
options_SA = optimoptions(@simulannealbnd, 'Display', 'iter');
lb = [-10, -10, -10];
ub = [10, 10, 10];
initial_guess = [1, 1, 1];

[params_optimized_SA, fval_SA] = simulannealbnd(objective_function, initial_guess, lb, ub, options_SA);

%Global optimisation using Genetic Algorithm
options_GA = optimoptions(@ga, 'Display', 'iter');
[params_optimized_GA, fval_GA] = ga(objective_function, length(initial_guess), [], [], [], [], lb, ub, [], options_GA);

%Global optimisation using Pattern Search
options_PS = optimoptions(@patternsearch, 'Display', 'iter');
[params_optimized_PS, fval_PS] = patternsearch(objective_function, initial_guess, [], [], [], [], lb, ub, [], options_PS);

%plotting the optimised models
y_optimized_SA = params_optimized_SA(1) * x.^2 + params_optimized_SA(2) * x + params_optimized_SA(3);
y_optimized_GA = params_optimized_GA(1) * x.^2 + params_optimized_GA(2) * x + params_optimized_GA(3);
y_optimized_PS = params_optimized_PS(1) * x.^2 + params_optimized_PS(2) * x + params_optimized_PS(3);

plot(x, y_optimized_SA, 'g', 'LineWidth', 2);
plot(x, y_optimized_GA, 'm', 'LineWidth', 2);
plot(x, y_optimized_PS, 'c', 'LineWidth', 2);

legend('Observed Data', 'True Model', 'SA', 'GA', 'PS');
hold off;

% Display optimized parameters and costs
disp('Optimized Parameters (Simulated Annealing):');
disp(params_optimized_SA);
disp(['Optimized Cost (Simulated Annealing): ', num2str(fval_SA)]);

disp('Optimized Parameters (Genetic Algorithm):');
disp(params_optimized_GA);
disp(['Optimized Cost (Genetic Algorithm): ', num2str(fval_GA)]);

disp('Optimized Parameters (Pattern Search):');
disp(params_optimized_PS);
disp(['Optimized Cost (Pattern Search): ', num2str(fval_PS)]);