rng(42);

global optimization_history;
optimization_history = [];

% Set paths to your plasma data files
csf_data_file1 = 'data/blattner2020_csf_concentration.csv';

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
    A = 84.523;
    r_p = 0.298;
    r_bp = 0.034;
    r_bc = 0.019;
    r_cp = 0.0154;
    sigma_A = 0.633;
    sigma_bp = 1.816;
    
    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
    	sigma_bc = params(1);
    	sigma_cp = params(2);
    	sigma_p = params(3);
    else
    	sigma_bc = 1.66;
        sigma_cp = 5.74;
        sigma_p = 3.61;
    end
    
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

function [c, ceq] = nonlinear_constraints(params)
    % Extract parameters
    r_bc = 0.019;
    r_bp = 0.034;
    sigma_bc = params(1);
    sigma_bp = 1.816;
    
    % Inequality constraints (c <= 0)
    % Wake constraint: r_bc + r_bp < 0.25 → r_bc + r_bp - 0.25 <= 0
    c(1) = r_bc + r_bp - 0.25;
    
    % Sleep constraint: sigma_bc*r_bc + sigma_bp*r_bp < 0.5 → sigma_bc*r_bc + sigma_bp*r_bp - 0.5 <= 0
    c(2) = sigma_bc * r_bc + sigma_bp * r_bp - 0.5;
    
    % No equality constraints
    ceq = [];
end

% Defining Objective Function for fmincon
% Updated objective function with steady-state penalty
function [total_error] = objective_function_fmincon(params, exp_csf1, fitting_error0)
    % Run simulation
    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [800,600,20], 0.01);
    
    % Extract last 36 hours of data
    brain_last_36hours_data = sol(233600:237200,1);  % Brain compartment (y(1))
    csf_last_36hours_data = sol(233600:237200,2);

    % Select data points corresponding to experimental time points
    time_indices = 1:200:3601;
    csf_model = csf_last_36hours_data(time_indices);

    % Ensure column vectors
    csf_model1 = csf_model(:);
    exp_csf1 = exp_csf1(:);

    % Calculate original fitting errors
    errors_C1 = (csf_model1 - exp_csf1).^ 2;

    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));

    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));

    fitting_error = nrmse_C1;

    % === UPDATED: Smooth, continuous steady-state penalty ===
    
    % Extract parameters
    r_bc = 0.019;
    r_bp = 0.034;
    sigma_bc = params(1);
    sigma_bp = 1.816;
    A = 84.523;
    sigma_A = 0.633;
    
    % Calculate theoretical steady-state values
    steady_state_wake = A / (r_bc + r_bp);
    steady_state_sleep = (sigma_A * A) / (sigma_bc * r_bc + sigma_bp * r_bp);
    
    % Create time vector for the analysis window
    start_time = (233600-1) * 0.01;
    time_analysis = start_time + (0:(length(brain_last_36hours_data)-1)) * 0.01;
    
    % Determine wake/sleep periods (wake: 8-24, sleep: 0-8)
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
    
    % Extract wake and sleep values
    wake_vals = brain_last_36hours_data(wake_mask);
    sleep_vals = brain_last_36hours_data(sleep_mask);
    
    % Avoid divide-by-zero: ensure steady_state > tiny
    eps_small = 1e-12;
    ss_w = max(abs(steady_state_wake), eps_small);
    ss_s = max(abs(steady_state_sleep), eps_small);
    
    % Calculate relative deviations from steady state
    wake_dev = abs(wake_vals - steady_state_wake) / ss_w;
    sleep_dev = abs(sleep_vals - steady_state_sleep) / ss_s;
    
    % Continuous closeness penalty parameters
    penalty_threshold = 0.10;  % 10% tolerance (recommended default)
    
    % Calculate continuous closeness scores (0 when dev >= threshold, 1 when dev = 0)
    cw = max(0, (penalty_threshold - wake_dev) / penalty_threshold);
    cs = max(0, (penalty_threshold - sleep_dev) / penalty_threshold);
    
    % Raw penalty (mean squared closeness), normalized between 0 and 1
    penalty_raw_w = mean(cw.^2);
    penalty_raw_s = mean(cs.^2);
    penalty_raw = 0.5 * (penalty_raw_w + penalty_raw_s);
    
    % Auto-calibrate penalty weight so penalty contributes ~5% initially
    p_frac = 0.05;  % aim for penalty to be 5% of fitting_error initially
    
    % Use provided fitting_error0 or fallback to current fitting_error
    if nargin < 3 || isempty(fitting_error0)
        fitting_error0 = fitting_error + 1e-12;
    end
    
    if penalty_raw > 0
        penalty_weight = p_frac * fitting_error0 / penalty_raw;
    else
        penalty_weight = 1e-6;
    end
    
    % Cap weight to avoid absurdly large values (optional)
    max_weight = 100;
    penalty_weight = min(penalty_weight, max_weight);
    
    % Final steady-state penalty
    steady_state_penalty = penalty_weight * penalty_raw;
    
    % Combined objective function
    total_error = fitting_error + steady_state_penalty;

    % Print current parameters and errors for debugging
    fprintf('Parameters: s_bc=%.5f, s_cp=%.5f, s_p=%.5f\n', ...
            params(1), params(2), params(3));
    fprintf('Individual Study Errors: CSF1=%.5F\n', ...
            error_C1);
    fprintf('Fitting Error: %.5f, Steady-State Penalty: %.5f (weight=%.5f), Total: %.5f\n', ...
            fitting_error, steady_state_penalty, penalty_weight, total_error);
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
    params_matrix = zeros(num_starts, 3);
    fvals = zeros(num_starts, 1);
    for i = 1:num_starts
        params_matrix(i,:) = solutions(i).params;
        fvals(i) = solutions(i).Fval;
    end
    
    % Create parameter names for plotting
    param_names = {'s_bc', 's_cp', 's_p'};
    
    % Get colormap with num_starts distinct colors
    colors = jet(num_starts);
    
    % Create figure with subplots for each parameter
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot parameters in subplots
    for i = 1:3
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
    saveas(gcf, 'blattner_fit/parameter_space.png');
    print('blattner_fit/parameter_space_highres', '-dpng', '-r300');
    
    % NEW CODE: Save parameter space data to CSV files
    for param_idx = 1:3
        % Create and open CSV file for this parameter
        param_name = param_names{param_idx};
        csv_filename = sprintf('blattner_fit/parameter_space_%s.csv', param_name);
        
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
lb = [0, 0, 0];
ub = [1.66, 5.74, 3.61];
initial_guess = 0.5 * (lb + ub);

fitting_error0 = objective_function_fmincon(initial_guess, csf_conc_exp1);
fprintf('Initial fitting error for penalty calibration: %.6f\n', fitting_error0);

% Create figure for visualization
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

% Set up the optimization problem with fmincon
problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function_fmincon(params, csf_conc_exp1), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'nonlcon', @nonlinear_constraints, ...
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
if ~exist('blattner_fit', 'dir')
    mkdir('blattner_fit');
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
    save(sprintf('blattner_fit/run_%d_results.mat', i), ...
        'run_results');
    
    % Log the results
    fprintf('\nRun %d:\n', i);
    fprintf('Parameters: s_bc=%.5f, s_cp=%.5f, s_p=%.5f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3));
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
param_names = {'s_bc', 's_cp', 's_p'};
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
fileID_corr = fopen('blattner_fit/correlation_matrix.txt', 'w');
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
save('blattner_fit/final_results.mat', 'final_results');

% Save best parameters with full precision in a MAT file
best_params_struct = struct();
param_names = {'s_bc', 's_cp', 's_p'};
for i = 1:length(param_names)
    best_params_struct.(param_names{i}) = best_params(i);
end
save('blattner_fit/best_parameters.mat', 'best_params_struct');

% Also save just the raw best_params array if needed
save('blattner_fit/best_params_array.mat', 'best_params');

% Create a text file to save results
fileID = fopen('blattner_fit/detailed_results.txt', 'w');

% Write header
fprintf(fileID, 'Optimization Results Summary\n');
fprintf(fileID, '===========================\n\n');

% For each run
for i = 1:num_starts
    fprintf(fileID, 'Run %d:\n', i);
    fprintf('Parameters: s_bc=%.5f, s_cp=%.5f, s_p=%.5f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3));
    
    % Calculate individual errors for this solution
    [t_run, sol_run] = euler(@(t,y) model(t,y,solutions(i).params), [0, 24*100], [800,600,20], 0.01);
    
    % Extract last 36 hours and calculate errors (similar to your existing error calculation)
    brain_last_36hours_data = sol_run(233600:237200,1);
    csf_last_36hours = sol_run(233600:237200,2);
    
    time_indices = 1:200:3601;
    csf_model = csf_last_36hours(time_indices);
    
    % Ensure column vectors
    csf_model1 = csf_model(:);
    exp_csf1 = csf_conc_exp1(:);

    errors_C1 = (csf_model1 - exp_csf1).^ 2;

    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));
    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));

    fitting_error = nrmse_C1;
    
    % === NEW: Calculate steady-state penalty ===
    
    % Extract parameters
    r_bc = 0.019;
    r_bp = 0.034;
    sigma_bc = solutions(i).params(1);
    sigma_bp = 1.816;
    A = 84.523;
    sigma_A = 0.633;
    
    % Calculate theoretical steady-state values
    steady_state_wake = A / (r_bc + r_bp);
    steady_state_sleep = (sigma_A * A) / (sigma_bc * r_bc + sigma_bp * r_bp);
    
    % Create time vector for the analysis window using your existing indexing
    % Indices 233600:237600 correspond to times (233600-1)*0.01 to (237600-1)*0.01
    start_time = (233600-1) * 0.01;  % = 2335.99 hours
    time_analysis = start_time + (0:(length(brain_last_36hours_data)-1)) * 0.01;
    
    % Determine wake/sleep periods (wake: 8-24, sleep: 0-8) for each time point
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
    
    % Calculate penalties for being too close to steady state
    penalty_threshold = 0.05; % 5% tolerance - adjust as needed
    penalty_weight = 0.5;     % Weight of penalty relative to fitting error - adjust as needed
    
    % Wake period penalty
    wake_brain_values = brain_last_36hours_data(wake_mask);
    wake_deviations = abs(wake_brain_values - steady_state_wake) / steady_state_wake;
    wake_penalty_mask = wake_deviations < penalty_threshold;
    wake_penalty = sum(wake_penalty_mask) / length(wake_brain_values);
    
    % Sleep period penalty
    sleep_brain_values = brain_last_36hours_data(sleep_mask);
    sleep_deviations = abs(sleep_brain_values - steady_state_sleep) / steady_state_sleep;
    sleep_penalty_mask = sleep_deviations < penalty_threshold;
    sleep_penalty = sum(sleep_penalty_mask) / length(sleep_brain_values);
    
    % Total steady-state penalty
    steady_state_penalty = penalty_weight * (wake_penalty + sleep_penalty) / 2;
    
    % Combined objective function
    total_error = fitting_error + steady_state_penalty;

    % Write errors to file
    fprintf(fileID, 'Individual Errors:\n');
    fprintf(fileID, '  CSF1 Error: %.6f\n', error_C1);
    fprintf(fileID, 'Total Error: %.6f\n', total_error);
end

% Write best solution
fprintf(fileID, '\nBest Solution:\n');
fprintf(fileID, 'Parameters: s_bc=%.5f, s_cp=%.5f, s_p=%.5f\n', ...
    best_params(1), best_params(2), best_params(3));
fprintf(fileID, 'Objective value: %.6f\n', best_fval);

fclose(fileID);

% Save loss curves for each run
for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = loss_curve(:, end); % Last column has the fval
    loss_data = [iterations(:), loss_values(:)];
    writematrix(loss_data, sprintf('blattner_fit/loss_curve_run_%d.csv', i));
end

% Plot loss curves for each run
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
saveas(gcf, 'blattner_fit/loss_curves.png');

% Plot parameter space
plot_parameter_space(solutions, num_starts);

% Print final best results
fprintf('\nBest solution found:\n');
fprintf('Parameters: s_bc=%.5f, s_cp=%.5f, s_p=%.5f\n', ...
    best_params(1), best_params(2), best_params(3));
