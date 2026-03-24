rng(42);

global optimization_history;
optimization_history = [];

% Data paths
csf_data_file3 = 'data/liu_csf_wake_conc.csv';
plasma_data_file1 = 'data/liu_plasma_wake_conc.csv';

csf_data3 = readtable(csf_data_file3);
plasma_data1 = readtable(plasma_data_file1);

% Extract
time_exp3 = csf_data3.Time;
csf_conc_exp3 = csf_data3.Concentration;
plasma_conc_exp1 = plasma_data1.Concentration;
csf_lsd3 = csf_data3.LSD;
csf_usd3 = csf_data3.USD;
plasma_lsd1 = plasma_data1.LSD;
plasma_usd1 = plasma_data1.USD;
csf_std3 = (csf_usd3 - csf_lsd3)/2;
plasma_std1 = (plasma_usd1 - plasma_lsd1)/2;

% Model
function dydt_n = model(t, y, params)
    r_bc = params(1);
    r_bp = params(2);
    r_cp = params(3);
    sigma_bc = params(4);
    sigma_bp = params(5);
    sigma_cp = params(6);
    sigma_p = params(7);
    A = params(8);
    sigma_A = params(9);
    r_p = params(10);

    % Switch
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

% Obj func
function [total_error] = objective_function_fmincon(params, exp_csf1, exp_csf2, exp_csf3, exp_plasma1, fitting_error0)
    % Run simulation
    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [800,600,20], 0.01);
    
    brain_last_36hours_data = sol(233600:237600,1);
    csf_last_36hours_data = sol(233600:237600,2);
    plasma_last_36hours_data = sol(233600:237600,3);

    time_indices = 1:200:4001;
    csf_model = csf_last_36hours_data(time_indices);
    plasma_model = plasma_last_36hours_data(time_indices);

    selected_indices3 = [1:9, 13:21];
    csf_model3 = csf_model(selected_indices3);
    plasma_model1 = plasma_model(selected_indices3);

    csf_model3 = csf_model3(:);
    plasma_model1 = plasma_model1(:);
    exp_csf3 = exp_csf3(:);
    exp_plasma1 = exp_plasma1(:);

    errors_C3 = (csf_model3 - exp_csf3).^ 2;
    errors_P1 = (plasma_model1 - exp_plasma1).^2;

    error_C3 = sqrt(sum(errors_C3) / length(exp_csf3));
    error_P1 = sqrt(sum(errors_P1) / length(exp_plasma1));

    nrmse_C3 = error_C3/abs(max(exp_csf3) - min(exp_csf3));
    nrmse_P1 = error_P1/abs(max(exp_plasma1) - min(exp_plasma1));

    fitting_error = (nrmse_C3+nrmse_P1)/2;

    % Penalty
    % Extract parameters
    r_bc = params(1);
    r_bp = params(2);
    sigma_bc = params(4);
    sigma_bp = params(5);
    A = params(8);
    sigma_A = params(9);
    
    % Theoretical steady-state values
    steady_state_wake = A / (r_bc + r_bp);
    steady_state_sleep = (sigma_A * A) / (sigma_bc * r_bc + sigma_bp * r_bp);

    start_time = (233600-1) * 0.01;
    time_analysis = start_time + (0:(length(brain_last_36hours_data)-1)) * 0.01;

    wake_mask = false(size(brain_last_36hours_data));
    sleep_mask = false(size(brain_last_36hours_data));
    
    for i = 1:length(time_analysis)
        hour_of_day = mod(time_analysis(i), 24);
        if hour_of_day >= 8 && hour_of_day < 24
            wake_mask(i) = true;
        else
            sleep_mask(i) = true;
        end
    end
    
    wake_vals = brain_last_36hours_data(wake_mask);
    sleep_vals = brain_last_36hours_data(sleep_mask);
    eps_small = 1e-12;
    ss_w = max(abs(steady_state_wake), eps_small);
    ss_s = max(abs(steady_state_sleep), eps_small);

    wake_dev = abs(wake_vals - steady_state_wake) / ss_w;
    sleep_dev = abs(sleep_vals - steady_state_sleep) / ss_s;

    penalty_threshold = 0.10;  % 10% tolerance

    cw = max(0, (penalty_threshold - wake_dev) / penalty_threshold);
    cs = max(0, (penalty_threshold - sleep_dev) / penalty_threshold);

    penalty_raw_w = mean(cw.^2);
    penalty_raw_s = mean(cs.^2);
    penalty_raw = 0.5 * (penalty_raw_w + penalty_raw_s);

    p_frac = 0.05;  % penalty to be 5% of fitting_error

    if nargin < 6 || isempty(fitting_error0)
        fitting_error0 = fitting_error + 1e-12;
    end
    
    if penalty_raw > 0
        penalty_weight = p_frac * fitting_error0 / penalty_raw;
    else
        penalty_weight = 1e-6;
    end

    max_weight = 100;
    penalty_weight = min(penalty_weight, max_weight);

    steady_state_penalty = penalty_weight * penalty_raw;
    total_error = fitting_error + steady_state_penalty;

    fprintf('Parameters: r_bc=%.5f, r_bp=%.5f, r_cp=%.5f, sigma_bc=%.5f, sigma_bp=%.5f, sigma_cp=%.5f, sigma_p=%.5f, A=%.5f, sigma_A=%.5f, r_p=%.5f\n', ...
            params(1), params(2), params(3), params(4), params(5), params(6), params(7), params(8), params(9), params(10));
    fprintf('Individual Study Errors: CSF3=%.5f, Plasma=%.5f\n', ...
            error_C3, error_P1);
    fprintf('Fitting Error: %.5f, Steady-State Penalty: %.5f (weight=%.5f), Total: %.5f\n', ...
            fitting_error, steady_state_penalty, penalty_weight, total_error);
    fprintf('Wake SS: %.5f, Sleep SS: %.5f, Penalty Raw: %.5f\n', ...
            steady_state_wake, steady_state_sleep, penalty_raw);
end

function stop = output_function(x, optimValues, state)
    global optimization_history;
    stop = false;
    
    if strcmp(state, 'iter')
        optimization_history = [optimization_history; x];
    end
end

function plot_parameter_space(solutions, num_starts)
    params_matrix = zeros(num_starts, 10);
    fvals = zeros(num_starts, 1);
    for i = 1:num_starts
        params_matrix(i,:) = solutions(i).params;
        fvals(i) = solutions(i).Fval;
    end
    
    param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
    colors = jet(num_starts);
    figure('Position', [100, 100, 1200, 800]);
    for i = 1:10
        subplot(3, 4, i);
        hold on;

        for j = 1:num_starts
            scatter(params_matrix(j,i), fvals(j), 50, colors(j,:), 'filled');
        end
        
        xlabel(param_names{i});
        ylabel('Objective Value');
        title(sprintf('Effect of %s on Objective', param_names{i}));
        grid on;
        hold off;
    end
    
    subplot(3, 4, 11:12);
    hold on;

    h = zeros(num_starts, 1);
    for j = 1:num_starts
        h(j) = scatter(NaN, NaN, 50, colors(j,:), 'filled');
    end

    legend(h, arrayfun(@(x) sprintf('Run %d', x), 1:num_starts, 'UniformOutput', false), 'Location', 'bestoutside', 'NumColumns', ceil(num_starts/10));
    axis off;
    hold off;

    sgtitle('Parameter Space Exploration');
    saveas(gcf, 'liu_fit/parameter_space.png');
    print('liu_fit/parameter_space_highres', '-dpng', '-r300');

    for param_idx = 1:10
        param_name = param_names{param_idx};
        csv_filename = sprintf('liu_fit/parameter_space_%s.csv', param_name);
        
        csv_data = zeros(num_starts, 3);
        for run_idx = 1:num_starts
            csv_data(run_idx, 1) = run_idx;
            csv_data(run_idx, 2) = params_matrix(run_idx, param_idx);
            csv_data(run_idx, 3) = fvals(run_idx);
        end
        
        writematrix(csv_data, csv_filename);
        fprintf('Saved parameter space data for %s to %s\n', param_name, csv_filename);
    end
end

% ===============
% Main script
% ===============

% Bounds and initial guess from config
[lb, ub, initial_guess] = config();
num_starts = 20;

figure('Position', [100, 100, 1200, 800]);

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 2000, ...
    'Algorithm', 'interior-point', ...
    'FiniteDifferenceStepSize', 1e-4, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'OutputFcn', @output_function);

% Set problem
problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function_fmincon(params, csf_conc_exp3, plasma_conc_exp1), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'nonlcon', @nonlinear_constraints, ...
    'options', options);

% Custom start points
customStartPoints = zeros(num_starts, length(initial_guess));
customStartPoints(1,:) = initial_guess;

for i = 2:num_starts
    customStartPoints(i,:) = lb + rand(1, length(lb)) .* (ub - lb);
end

solutions = struct('params', cell(1, num_starts), 'Fval', cell(1, num_starts), 'Exitflag', cell(1, num_starts), 'loss_history', cell(1, num_starts));
all_trajectories = cell(num_starts, 1);
global optimization_history;

if ~exist('liu_fit', 'dir')
    mkdir('liu_fit');
end

best_fval = Inf;
best_params = [];
aic_values = zeros(num_starts, 1);

for i = 1:num_starts
    optimization_history = [];
    problem.x0 = customStartPoints(i,:);
    
    [x, fval, exitflag, output] = fmincon(problem);

    all_trajectories{i} = optimization_history;

    solutions(i).params = x(:)';
    solutions(i).Fval = fval;
    solutions(i).Exitflag = exitflag;

    if fval < best_fval
        best_fval = fval;
        best_params = x;
    end

    run_results = struct('params', x(:)', ...
        'fval', fval, ...
        'exitflag', exitflag, ...
        'trajectory', optimization_history);
    save(sprintf('liu_fit/run_%d_results.mat', i), ...
        'run_results');

    fprintf('\nRun %d:\n', i);
    fprintf('Parameters: r_bc=%.5f, r_bp=%.5f, r_cp=%.5f, sigma_bc=%.5f, sigma_bp=%.5f, sigma_cp=%.5f, sigma_p=%.5f, A=%.5f, sigma_A=%.5f, r_p=%.5f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3), solutions(i).params(4), solutions(i).params(5), solutions(i).params(6), solutions(i).params(7), ...
        solutions(i).params(8), solutions(i).params(9), solutions(i).params(10));
    fprintf('Objective value: %.6f\n', fval);
    fprintf('Completed run %d/%d\n', i, num_starts);
end

all_solutions = zeros(num_starts, length(initial_guess));
all_fvals = zeros(num_starts, 1);
for i = 1:num_starts
    all_solutions(i,:) = solutions(i).params;
    all_fvals(i) = solutions(i).Fval;
end

param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
correlation_matrix = corrcoef(all_solutions);

fprintf('\nParameter Correlation Matrix:\n');
fprintf('%-10s', 'Param');
for i = 1:length(param_names)
    fprintf('%-10s', param_names{i});
end
fprintf('\n');

for i = 1:length(param_names)
    fprintf('%-10s', param_names{i});
    for j = 1:length(param_names)
        fprintf('%-10.4f', correlation_matrix(i,j));
    end
    fprintf('\n');
end

fileID_corr = fopen('liu_fit/correlation_matrix.txt', 'w');
fprintf(fileID_corr, 'Parameter Correlation Matrix:\n');
fprintf(fileID_corr, '%-10s', 'Param');
for i = 1:length(param_names)
    fprintf(fileID_corr, '%-10s', param_names{i});
end
fprintf(fileID_corr, '\n');

for i = 1:length(param_names)
    fprintf(fileID_corr, '%-10s', param_names{i});
    for j = 1:length(param_names)
        fprintf(fileID_corr, '%-10.4f', correlation_matrix(i,j));
    end
    fprintf(fileID_corr, '\n');
end
fclose(fileID_corr);

final_results = struct('best_params', best_params, ...
    'best_fval', best_fval, ...
    'all_solutions', all_solutions, ...
    'all_fvals', all_fvals, ...
    'all_trajectories', all_trajectories, ...
    'correlation_matrix', correlation_matrix);
save('liu_fit/final_results.mat', 'final_results');

best_params_struct = struct();
param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
for i = 1:length(param_names)
    best_params_struct.(param_names{i}) = best_params(i);
end
save('liu_fit/best_parameters.mat', 'best_params_struct');
save('liu_fit/best_params_array.mat', 'best_params');
fileID = fopen('liu_fit/detailed_results.txt', 'w');

fprintf(fileID, 'Optimization Results Summary\n');
fprintf(fileID, '===========================\n\n');

% For each run
for i = 1:num_starts
    fprintf(fileID, 'Run %d:\n', i);
    fprintf(fileID, 'Parameters: r_bc=%.5f, r_bp=%.5f, r_cp=%.5f, sigma_bc=%.5f, sigma_bp=%.5f, sigma_cp=%.5f, sigma_p=%.5f, A=%.5f, sigma_A=%.5f, r_p=%.5f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3), solutions(i).params(4), solutions(i).params(5), solutions(i).params(6), solutions(i).params(7), ...
        solutions(i).params(8), solutions(i).params(9), solutions(i).params(10));
    
    [t_run, sol_run] = euler(@(t,y) model(t,y,solutions(i).params), [0, 24*100], [800,600,20], 0.01);

    brain_last_36hours_data = sol_run(233600:237600,1);
    csf_last_36hours = sol_run(233600:237600,2);
    plasma_last_36hours = sol_run(233600:237600,3);
    
    time_indices = 1:200:4001;
    csf_model = csf_last_36hours(time_indices);
    plasma_model = plasma_last_36hours(time_indices);

    selected_indices3 = [1:9, 13:21];
    csf_model3 = csf_model(selected_indices3);
    plasma_model1 = plasma_model(selected_indices3);

    csf_model3 = csf_model3(:);
    plasma_model1 = plasma_model1(:);
    exp_csf3 = csf_conc_exp3(:);
    exp_plasma1 = plasma_conc_exp1(:);

    errors_C3 = (csf_model3 - exp_csf3).^ 2;
    errors_P1 = (plasma_model1 - exp_plasma1).^2;

    error_C3 = sqrt(sum(errors_C3) / length(exp_csf3));
    error_P1 = sqrt(sum(errors_P1) / length(exp_plasma1));

    nrmse_C3 = error_C3/abs(max(exp_csf3) - min(exp_csf3));
    nrmse_P1 = error_P1/abs(max(exp_plasma1) - min(exp_plasma1));

    fitting_error = (nrmse_C3+nrmse_P1)/2;

    r_bc = solutions(i).params(1);
    r_bp = solutions(i).params(2);
    sigma_bc = solutions(i).params(4);
    sigma_bp = solutions(i).params(5);
    A = solutions(i).params(8);
    sigma_A = solutions(i).params(9);

    steady_state_wake = A / (r_bc + r_bp);
    steady_state_sleep = (sigma_A * A) / (sigma_bc * r_bc + sigma_bp * r_bp);
    
    start_time = (233600-1) * 0.01;
    time_analysis = start_time + (0:(length(brain_last_36hours_data)-1)) * 0.01;
    
    wake_mask = false(size(brain_last_36hours_data));
    sleep_mask = false(size(brain_last_36hours_data));
    
    for j = 1:length(time_analysis)
        hour_of_day = mod(time_analysis(j), 24);
        if hour_of_day >= 8 && hour_of_day < 24
            wake_mask(j) = true;
        else
            sleep_mask(j) = true;
        end
    end
    
    penalty_threshold = 0.05;
    penalty_weight = 0.5;
    
    wake_brain_values = brain_last_36hours_data(wake_mask);
    wake_deviations = abs(wake_brain_values - steady_state_wake) / steady_state_wake;
    wake_penalty_mask = wake_deviations < penalty_threshold;
    wake_penalty = sum(wake_penalty_mask) / length(wake_brain_values);
    
    sleep_brain_values = brain_last_36hours_data(sleep_mask);
    sleep_deviations = abs(sleep_brain_values - steady_state_sleep) / steady_state_sleep;
    sleep_penalty_mask = sleep_deviations < penalty_threshold;
    sleep_penalty = sum(sleep_penalty_mask) / length(sleep_brain_values);
    
    steady_state_penalty = penalty_weight * (wake_penalty + sleep_penalty) / 2;
    total_error = fitting_error + steady_state_penalty;

    fprintf(fileID, 'Individual Errors:\n');
    fprintf(fileID, '  CSF3 Error: %.6f\n', error_C3);
    fprintf(fileID, '  Plasma Error: %.6f\n', error_P1);
    fprintf(fileID, 'Total Error: %.6f\n', total_error);
end

% Write best solution
fprintf(fileID, '\nBest Solution:\n');
fprintf(fileID, 'Parameters: r_bc=%.5f, r_bp=%.5f, r_cp=%.5f, sigma_bc=%.5f, sigma_bp=%.5f, sigma_cp=%.5f, sigma_p=%.5f, A=%.5f, sigma_A=%.5f, r_p=%.5f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5), best_params(6), best_params(7), best_params(8), best_params(9), best_params(10));
fprintf(fileID, 'Objective value: %.6f\n', best_fval);

fclose(fileID);

for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = loss_curve(:, end);
    loss_data = [iterations(:), loss_values(:)];
    writematrix(loss_data, sprintf('liu_fit/loss_curve_run_%d.csv', i));
end

figure;
hold on;
for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = loss_curve(:, end);
    plot(iterations, loss_values, 'DisplayName', sprintf('Run %d', i));
end
xlabel('Iteration');
ylabel('Loss Value');
title('Loss Curves for Each Run');
legend;
hold off;
saveas(gcf, 'liu_fit/loss_curves.png');

plot_parameter_space(solutions, num_starts);

fprintf('\nBest solution found:\n');
fprintf('Parameters: r_bc=%.5f, r_bp=%.5f, r_cp=%.5f, sigma_bc=%.5f, sigma_bp=%.5f, sigma_cp=%.5f, sigma_p=%.5f, A=%.5f, sigma_A=%.5f, r_p=%.5f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5), best_params(6), best_params(7), best_params(8), best_params(9), best_params(10));
fprintf('Objective value: %.6f\n', best_fval);
