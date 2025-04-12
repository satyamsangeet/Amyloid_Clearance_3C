global optimization_history;
optimization_history = [];

% Set paths to your plasma data files
csf_data_file1 = 'data/blattner_wake_conc.csv';

% Read plasma data from both files
csf_data1 = readtable(csf_data_file1);

% Extract data from both plasma files
time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Concentration;
csf_lsd1 = csf_data1.LSD;
csf_usd1 = csf_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;

% Updated model function
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

% Defining Euler-Maruyama Method
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

% Defining Objective Function for fmincon
function [total_error] = objective_function_fmincon(params, exp_csf1)
    global optimization_history;
    optimization_history = [optimization_history; params];

    % Run simulation
    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
    
    % Extract last 36 hours of data
    csf_last_36hours_data = sol(233600:237400,2);

    % Select data points corresponding to experimental time points
    time_indices = 1:200:3801;
    csf_model = csf_last_36hours_data(time_indices);

    selected_indices1 = [1:8, 14:20];
    csf_model1 = csf_model(selected_indices1);
    
    % Ensure column vectors
    csf_model1 = csf_model1(:);
    exp_csf1 = exp_csf1(:);

    errors_C1 = (csf_model1 - exp_csf1) .^ 2;

    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));

    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));

    total_error = nrmse_C1;

    % Print current parameters and errors for debugging
    fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
            params(1), params(2), params(3), params(4), params(5), params(6), params(7), params(8), params(9), params(10));
    fprintf('Individual Study Errors: CSF1=%.4f\n', ...
            error_C1);
    fprintf('Total WRMSE: %.4f\n', total_error);
end

function stop = output_function(x, optimValues, state)
    global optimization_history;
    stop = false;
    
    if strcmp(state, 'iter')
        optimization_history = [optimization_history; x];
    end
end

% Function to plot parameter space
function plot_parameter_space(solutions, num_starts)
    % Extract parameters and objective values
    params_matrix = zeros(num_starts, 10);
    fvals = zeros(num_starts, 1);
    for i = 1:num_starts
        params_matrix(i,:) = solutions(i).params;
        fvals(i) = solutions(i).Fval;
    end
    
    % Create parameter names for plotting
    param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
    
    % Get colormap with num_starts distinct colors
    colors = jet(num_starts);
    
    % Create figure with subplots for each parameter
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot parameters in subplots
    for i = 1:10
        subplot(3, 4, i);
        hold on;
        
        % Plot each run's data point with a distinct color
        for j = 1:num_starts
            scatter(params_matrix(j,i), fvals(j), 50, colors(j,:), 'filled');
        end
        
        xlabel(param_names{i});
        ylabel('Objective Value');
        title(sprintf('Effect of %s on Objective', param_names{i}));
        grid on;
        hold off;
    end
    
    % Create a separate subplot for the legend
    subplot(3, 4, 11:12);
    hold on;
    
    % Create dummy points for legend
    h = zeros(num_starts, 1);
    for j = 1:num_starts
        h(j) = scatter(NaN, NaN, 50, colors(j,:), 'filled');
    end
    
    % Add legend with run numbers
    legend(h, arrayfun(@(x) sprintf('Run %d', x), 1:num_starts, 'UniformOutput', false), 'Location', 'bestoutside', 'NumColumns', ceil(num_starts/10));
    axis off;
    hold off;
    
    % Add title to the figure
    sgtitle('Parameter Space Exploration');
    
    % Save the figure
    saveas(gcf, 'global_all_params/parameter_space.png');
    print('global_all_params/parameter_space_highres', '-dpng', '-r300');
    
    % NEW CODE: Save parameter space data to CSV files
    for param_idx = 1:10
        % Create and open CSV file for this parameter
        param_name = param_names{param_idx};
        csv_filename = sprintf('global_all_params/parameter_space_%s.csv', param_name);
        
        % Prepare data in the required format: [run_number, parameter_value, objective_value]
        csv_data = zeros(num_starts, 3);
        for run_idx = 1:num_starts
            csv_data(run_idx, 1) = run_idx;  % Run number
            csv_data(run_idx, 2) = params_matrix(run_idx, param_idx);  % Parameter value
            csv_data(run_idx, 3) = fvals(run_idx);  % Objective function value
        end
        
        % Save to CSV
        writematrix(csv_data, csv_filename);
        fprintf('Saved parameter space data for %s to %s\n', param_name, csv_filename);
    end
end

% Main script
% Set optimization parameters
num_starts = 20;
initial_guess = [1.5, 0.3, 0.5, 2.5, 3, 3, 2, 11, 0.8, 0.28];
lb = [0.75, 0.01, 0.001, 1, 1, 1, 1, 9, 0.7, 0.23];
ub = [2.5, 1, 3, 4, 7, 7, 7, 14, 0.9, 0.34];

% Create figure for visualization
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

% Set up the optimization problem with fmincon
problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function_fmincon(params, csf_conc_exp1), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

% Generate custom start points
customStartPoints = zeros(num_starts, length(initial_guess));
customStartPoints(1,:) = initial_guess;

for i = 2:num_starts
    customStartPoints(i,:) = lb + rand(1, length(lb)) .* (ub - lb);
end

% Initialize arrays to store results with a regular struct array
solutions = struct('params', cell(1, num_starts), 'Fval', cell(1, num_starts), 'Exitflag', cell(1, num_starts), 'loss_history', cell(1, num_starts));
all_trajectories = cell(num_starts, 1);
global optimization_history;

% Create directory for results if it doesn't exist
if ~exist('global_all_params', 'dir')
    mkdir('global_all_params');
end

% Best solution tracking
best_fval = Inf;
best_params = [];
aic_values = zeros(num_starts, 1); % To store AIC values

% For each start point
for i = 1:num_starts
    % Reset optimization history for this run
    optimization_history = [];
    
    % Create problem with current start point
    problem.x0 = customStartPoints(i,:);
    
    % Run optimization
    [x, fval, exitflag, output] = fmincon(problem);
    
    % Store trajectory
    all_trajectories{i} = optimization_history;
    
    % Store solution in struct
    solutions(i).params = x(:)';  % Store optimized parameters in 'params'
    solutions(i).Fval = fval;
    solutions(i).Exitflag = exitflag;
    
    % Update best solution if current is better
    if fval < best_fval
        best_fval = fval;
        best_params = x;
    end
    
    % Save individual run results
    run_results = struct('params', x(:)', ...
        'fval', fval, ...
        'exitflag', exitflag, ...
        'trajectory', optimization_history);
    save(sprintf('global_all_params/run_%d_results.mat', i), ...
        'run_results');
    
    % Log the results
    fprintf('\nRun %d:\n', i);
    fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3), solutions(i).params(4), solutions(i).params(5), solutions(i).params(6), solutions(i).params(7), ...
        solutions(i).params(8), solutions(i).params(9), solutions(i).params(10));
    fprintf('Objective value: %.6f\n', fval);
    fprintf('Completed run %d/%d\n', i, num_starts);
end

% Extract all solutions and function values
all_solutions = zeros(num_starts, length(initial_guess));
all_fvals = zeros(num_starts, 1);
for i = 1:num_starts
    all_solutions(i,:) = solutions(i).params;
    all_fvals(i) = solutions(i).Fval;
end

% Calculate correlation matrix for parameters
param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
correlation_matrix = corrcoef(all_solutions);

% Print correlation matrix
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

% Also save correlation matrix to file
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

% Save overall results
final_results = struct('best_params', best_params, ...
    'best_fval', best_fval, ...
    'all_solutions', all_solutions, ...
    'all_fvals', all_fvals, ...
    'all_trajectories', all_trajectories, ...
    'correlation_matrix', correlation_matrix);
save('global_all_params/final_results.mat', 'final_results');

% Save best parameters with full precision in a MAT file
best_params_struct = struct();
param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
for i = 1:length(param_names)
    best_params_struct.(param_names{i}) = best_params(i);
end
save('global_all_params/best_parameters.mat', 'best_params_struct');

% Also save just the raw best_params array if needed
save('global_all_params/best_params_array.mat', 'best_params');

% Create a text file to save results
fileID = fopen('global_all_params/detailed_results.txt', 'w');

% Write header
fprintf(fileID, 'Optimization Results Summary\n');
fprintf(fileID, '===========================\n\n');

% For each run
for i = 1:num_starts
    fprintf(fileID, 'Run %d:\n', i);
    fprintf(fileID, 'Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3), solutions(i).params(4), solutions(i).params(5), solutions(i).params(6), solutions(i).params(7), ...
        solutions(i).params(8), solutions(i).params(9), solutions(i).params(10));
    
    % Calculate individual errors for this solution
    [t_run, sol_run] = euler(@(t,y) model(t,y,solutions(i).params), [0, 24*100], [0,600,15.5], 0.01);
    
    % Extract last 36 hours and calculate errors (similar to your existing error calculation)
    csf_last_36hours = sol_run(233600:237400,2);
    
    time_indices = 1:200:3801;
    csf_model = csf_last_36hours(time_indices);
    
    selected_indices1 = [1:8, 14:20];
    csf_model1 = csf_model(selected_indices1);
    
    % Ensure column vectors
    csf_model1 = csf_model1(:);
    exp_csf1 = csf_conc_exp1(:);

    errors_C1 = (csf_model1 - exp_csf1) .^ 2;

    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));

    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));

    total_error = nrmse_C1;

    % Write errors to file
    fprintf(fileID, 'Individual Errors:\n');
    fprintf(fileID, '  CSF1 Error: %.6f\n', error_C1);
    fprintf(fileID, 'Total Error: %.6f\n', total_error);
end

% Write best solution
fprintf(fileID, '\nBest Solution:\n');
fprintf(fileID, 'Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5), best_params(6), best_params(7), best_params(8), best_params(9), best_params(10));
fprintf(fileID, 'Objective value: %.6f\n', best_fval);

fclose(fileID);

% Save loss curves for each run
for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = arrayfun(@(idx) objective_function_fmincon(loss_curve(idx,:), csf_conc_exp1), iterations);
    loss_data = [iterations(:), loss_values(:)];
    writematrix(loss_data, sprintf('global_all_params/loss_curve_run_%d.csv', i));
end

% Plot loss curves for each run
figure;
hold on;
for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = arrayfun(@(idx) objective_function_fmincon(loss_curve(idx,:), csf_conc_exp1), iterations);
    plot(iterations, loss_values, 'DisplayName', sprintf('Run %d', i));
end
xlabel('Iteration');
ylabel('Loss Value');
title('Loss Curves for Each Run');
legend;
hold off;
saveas(gcf, 'global_all_params/loss_curves.png');

% Plot parameter space
plot_parameter_space(solutions, num_starts);

% Print final best results
fprintf('\nBest solution found:\n');
fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4), best_params(5), best_params(6), best_params(7), best_params(8), best_params(9), best_params(10));
fprintf('Objective value: %.6f\n', best_fval);