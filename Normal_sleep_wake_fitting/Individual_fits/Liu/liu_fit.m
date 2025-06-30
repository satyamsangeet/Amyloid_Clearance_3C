rng(42);
global optimization_history;
optimization_history = [];

csf_data_file1 = 'data/liu_csf_wake_conc.csv';
plasma_data_file1 = 'data/liu_plasma_wake_conc.csv';
csf_data1 = readtable(csf_data_file1);
plasma_data1 = readtable(plasma_data_file1);
time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Concentration;
plasma_conc_exp1 = plasma_data1.Concentration;
csf_lsd1 = csf_data1.LSD;
csf_usd1 = csf_data1.USD;
plasma_lsd1 = plasma_data1.LSD;
plasma_usd1 = plasma_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;
plasma_std1 = (plasma_usd1 - plasma_lsd1)/2;

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

    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);

    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function [t, w] = euler(F, endpoints, initial_conditions, ts)
    if length(endpoints) == 2
        h = ts; %delta_t (seconds)
        total_time = endpoints(2) - endpoints(1);
        num_steps = floor(total_time / h);
        t = linspace(endpoints(1), endpoints(2), num_steps + 1); %Creating time vector
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints;
    end
    w = zeros(num_steps+1, length(initial_conditions));
    w(1,:) = initial_conditions;
    for k = 1:num_steps
        w(k+1,:) = w(k,:) + F(t(k), w(k,:))' * h;
    end
    t = t(:);
end

function [total_error] = objective_function_fmincon(params, exp_csf1, exp_plasma1)
    global optimization_history;
    optimization_history = [optimization_history; params];

    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
    
    csf_last_36hours_data = sol(233600:237400,2);
    plasma_last_36hours_data = sol(233600:237400,3);
    time_indices = 1:200:3801;
    csf_model = csf_last_36hours_data(time_indices);
    plasma_model = plasma_last_36hours_data(time_indices);

    selected_indices1 = [1:7, 14:20];
    csf_model1 = csf_model(selected_indices1);
    plasma_model1 = plasma_model(selected_indices1);

    csf_model1 = csf_model1(:);
    plasma_model1 = plasma_model1(:);
    exp_csf1 = exp_csf1(:);
    exp_plasma1 = exp_plasma1(:);

    errors_C1 = (csf_model1 - exp_csf1).^ 2;
    errors_P1 = (plasma_model1 - exp_plasma1).^2;
    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));
    error_P1 = sqrt(sum(errors_P1) / length(exp_plasma1));
    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));
    nrmse_P1 = error_P1/abs(max(exp_plasma1) - min(exp_plasma1));
    total_error = (nrmse_C1+nrmse_P1)/2;

    fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
            params(1), params(2), params(3), params(4), params(5), params(6), params(7), params(8), params(9), params(10));
    fprintf('Individual Study Errors: CSF1=%.4f, Plasma=%.4f\n', error_C1, error_P1);
    fprintf('Total WRMSE: %.4f\n', total_error);
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
    
    %saveas(gcf, 'global_all_params/parameter_space.png');
    %print('global_all_params/parameter_space_highres', '-dpng', '-r300');
    
    for param_idx = 1:10
        param_name = param_names{param_idx};
        csv_filename = sprintf('global_all_params/parameter_space_%s.csv', param_name);
        csv_data = zeros(num_starts, 3);
        for run_idx = 1:num_starts
            csv_data(run_idx, 1) = run_idx;  % Run number
            csv_data(run_idx, 2) = params_matrix(run_idx, param_idx);  % Parameter value
            csv_data(run_idx, 3) = fvals(run_idx);  % Objective function value
        end

        writematrix(csv_data, csv_filename);
        fprintf('Saved parameter space data for %s to %s\n', param_name, csv_filename);
    end
end

num_starts = 20;
initial_guess = [1.5, 0.3, 0.5, 2.5, 3, 3, 2, 11, 0.8, 0.28];
lb = [0.75, 0.01, 0.001, 1, 1, 1, 1, 9, 0.7, 0.23];
ub = [2.5, 1, 3, 4, 7, 7, 7, 14, 0.9, 0.34];

figure('Position', [100, 100, 1200, 800]);
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'Algorithm', 'interior-point', ...
    'FiniteDifferenceStepSize', 1e-4, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'OutputFcn', @output_function);

problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function_fmincon(params, csf_conc_exp1, plasma_conc_exp1), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

customStartPoints = zeros(num_starts, length(initial_guess));
customStartPoints(1,:) = initial_guess;

for i = 2:num_starts
    customStartPoints(i,:) = lb + rand(1, length(lb)) .* (ub - lb);
end

solutions = struct('params', cell(1, num_starts), 'Fval', cell(1, num_starts), 'Exitflag', cell(1, num_starts), 'loss_history', cell(1, num_starts));
all_trajectories = cell(num_starts, 1);
global optimization_history;

if ~exist('global_all_params', 'dir')
    mkdir('global_all_params');
end

best_fval = Inf;
best_params = [];

for i = 1:num_starts
    optimization_history = [];
    problem.x0 = customStartPoints(i,:);
    [x, fval, exitflag, output] = fmincon(problem);

    all_trajectories{i} = optimization_history;
    solutions(i).params = x(:)';  % Store optimized parameters in 'params'
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
    save(sprintf('global_all_params/run_%d_results.mat', i), ...
        'run_results');
    
    fprintf('\nRun %d:\n', i);
    fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
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

fileID_corr = fopen('global_all_params/correlation_matrix.txt', 'w');
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
save('global_all_params/final_results.mat', 'final_results');

best_params_struct = struct();
param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
for i = 1:length(param_names)
    best_params_struct.(param_names{i}) = best_params(i);
end
save('global_all_params/best_parameters.mat', 'best_params_struct');
save('global_all_params/best_params_array.mat', 'best_params');

fileID = fopen('global_all_params/detailed_results.txt', 'w');
fprintf(fileID, 'Optimization Results Summary\n');
fprintf(fileID, '===========================\n\n');

for i = 1:num_starts
    fprintf(fileID, 'Run %d:\n', i);
    fprintf(fileID, 'Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3), solutions(i).params(4), solutions(i).params(5), solutions(i).params(6), solutions(i).params(7), ...
        solutions(i).params(8), solutions(i).params(9), solutions(i).params(10));
    
    [t_run, sol_run] = euler(@(t,y) model(t,y,solutions(i).params), [0, 24*100], [0,600,15.5], 0.01);
    
    csf_last_36hours = sol_run(233600:237400,2);
    plasma_last_36hours = sol_run(233600:237400,3);
    
    time_indices = 1:200:3801;
    csf_model = csf_last_36hours(time_indices);
    plasma_model = plasma_last_36hours(time_indices);
    
    selected_indices1 = [1:7, 14:20];
    csf_model1 = csf_model(selected_indices1);
    plasma_model1 = plasma_model(selected_indices1);

    csf_model1 = csf_model1(:);
    plasma_model1 = plasma_model1(:);
    exp_csf1 = csf_conc_exp1(:);
    exp_plasma1 = plasma_conc_exp1(:);

    errors_C1 = (csf_model1 - exp_csf1).^ 2;
    errors_P1 = (plasma_model1 - exp_plasma1).^2;
    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));
    error_P1 = sqrt(sum(errors_P1) / length(exp_plasma1));
    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));
    nrmse_P1 = error_P1/abs(max(exp_plasma1) - min(exp_plasma1));
    total_error = (nrmse_C1+nrmse_P1)/2;

    fprintf(fileID, 'Individual Errors:\n');
    fprintf(fileID, '  CSF1 Error: %.6f\n', error_C1);
    fprintf(fileID, '  Plasma Error: %.6f\n', error_P1);
    fprintf(fileID, 'Total Error: %.6f\n', total_error);
end

fprintf(fileID, '\nBest Solution:\n');
fprintf(fileID, 'Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5), best_params(6), best_params(7), best_params(8), best_params(9), best_params(10));
fprintf(fileID, 'Objective value: %.6f\n', best_fval);
fclose(fileID);

for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = arrayfun(@(idx) objective_function_fmincon(loss_curve(idx,:), csf_conc_exp1, plasma_conc_exp1), iterations);
    loss_data = [iterations(:), loss_values(:)];
    writematrix(loss_data, sprintf('global_all_params/loss_curve_run_%d.csv', i));
end

figure;
hold on;
for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = arrayfun(@(idx) objective_function_fmincon(loss_curve(idx,:), csf_conc_exp1, plasma_conc_exp1), iterations);
    plot(iterations, loss_values, 'DisplayName', sprintf('Run %d', i));
end
xlabel('Iteration');
ylabel('Loss Value');
title('Loss Curves for Each Run');
legend;
hold off;
%saveas(gcf, 'global_all_params/loss_curves.png');

plot_parameter_space(solutions, num_starts);

fprintf('\nBest solution found:\n');
fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5), best_params(6), best_params(7), best_params(8), best_params(9), best_params(10));
fprintf('Objective value: %.6f\n', best_fval);
