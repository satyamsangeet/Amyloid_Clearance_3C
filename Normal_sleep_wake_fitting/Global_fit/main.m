clear;
clc;

rng(42);

global optimization_history;

% Loading data
[data, time] = load_data();

% Optimization requirements
num_starts = 20;
initial_guess = [1.5, 0.3, 0.5, 2.5, 3, 3, 2, 11, 0.8, 0.28];
lb = [0.75, 0.01, 0.001, 1, 1, 1, 1, 9, 0.7, 0.23];
ub = [2.5, 1, 3, 4, 7, 7, 7, 14, 0.9, 0.34];

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'Algorithm', 'interior-point', ...
    'FiniteDifferenceStepSize', 1e-4, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'OutputFcn', @output_function);

% Creating optimization problem
problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function(params, data.csf_conc_exp1, data.csf_conc_exp2, data.csf_conc_exp3, data.plasma_conc_exp1), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

% Generating the custom start points
customStartPoints = generate_start_points(initial_guess, lb, ub, num_starts);

% Run optim, save
[solutions, best_params, best_fval, all_trajectories] = run_optimization(problem, customStartPoints, num_starts, data);

plot_parameter_space(solutions, num_starts);

% Print
fprintf('\nBest solution found:\n');
fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5), best_params(6), best_params(7), best_params(8), best_params(9), best_params(10));
fprintf('Objective value: %.6f\n', best_fval);
