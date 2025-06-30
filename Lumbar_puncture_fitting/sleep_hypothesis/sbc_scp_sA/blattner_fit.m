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
    A = 13.791866842295299;
    r_p = 0.285014821392297;
    r_bc = 2.349030412757668;
    r_bp = 0.092859916742039;
    r_cp = 0.006556023403031;
    sigma_p = 5.859219697739381;
    sigma_bp = 2.526995585282243;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        sigma_bc = params(1);
        sigma_cp = params(2);
        sigma_A = params(3);
    else
        sigma_A = 0.718493031801368;
    	sigma_bc = 1.117763182229140;
    	sigma_cp = 6.416831718871831;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    % Define the ODE system
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
    csf_last_36hours_data = sol(233600:237200,2);

    % Select data points corresponding to experimental time points
    time_indices = 1:200:3601;
    csf_model = csf_last_36hours_data(time_indices);
    
    % Ensure column vectors
    csf_model1 = csf_model(:);
    exp_csf1 = exp_csf1(:);

    errors_C1 = (csf_model1 - exp_csf1) .^ 2;
    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));
    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));

    total_error = nrmse_C1;

    % Print current parameters and errors for debugging
    fprintf('Parameters: s_bc=%.5f, s_cp=%.5f, s_A=%.5f\n', ...
            params(1), params(2), params(3));
    fprintf('Individual Study Errors: CSF1=%.5f\n', ...
            error_C1);
    fprintf('Total WRMSE: %.5f\n', total_error);
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
    param_names = {'s_bc', 's_cp', 's_A'};
    
    % Get colormap with num_starts distinct colors
    colors = jet(num_starts);
    
    % Create figure with subplots for each parameter
    figure('Position', [100, 100, 1200, 800]);
    
    % Create a cell array to store plot handles for legend
    plot_handles = cell(num_starts, 1);
    
    for i = 1:3
        subplot(3, 4, i);
        hold on;
        
        % Plot each run's data point with a distinct color
        for j = 1:num_starts
            h = scatter(params_matrix(j,i), fvals(j), 50, colors(j,:), 'filled');
            
            % Store handle for first subplot to use in legend
            if i == 1
                plot_handles{j} = h;
            end
        end
        
        xlabel(param_names{i});
        ylabel('Objective Value');
        title(sprintf('Effect of %s on Objective', param_names{i}));
        grid on;
        hold off;
    end
    
    % Add legend to the last subplot
    subplot(3, 4, 11);  % Use position 11 for legend
    hold on;
    
    % Create dummy invisible plot for legend
    for j = 1:num_starts
        scatter(NaN, NaN, 50, colors(j,:), 'filled', 'DisplayName', sprintf('Run %d', j));
    end
    
    legend('Location', 'best');
    axis off;
    sgtitle('Parameter Space Exploration');
    saveas(gcf, 'blattner_exp_fit3/parameter_space.png');
    print('blattner_exp_fit3/parameter_space_highres', '-dpng', '-r300');

    % NEW CODE: Save parameter space data to CSV files
    for param_idx = 1:3
        % Create and open CSV file for this parameter
        param_name = param_names{param_idx};
        csv_filename = sprintf('blattner_exp_fit3/parameter_space_%s.csv', param_name);
        
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
initial_guess = [0.8, 5, 0.8];
lb = [0.01, 0.01, 0.718];
ub = [1.117, 6.416, 0.9];

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
if ~exist('blattner_exp_fit3', 'dir')
    mkdir('blattner_exp_fit3');
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
    save(sprintf('blattner_exp_fit3/run_%d_results.mat', i), ...
        'run_results');
    
    % Log the results
    fprintf('\nRun %d:\n', i);
    fprintf('Parameters: s_bc=%.5f, s_cp=%.5f, s_A=%.5f\n', ...
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

% Save overall results
final_results = struct('best_params', best_params, ...
    'best_fval', best_fval, ...
    'all_solutions', all_solutions, ...
    'all_fvals', all_fvals, ...
    'all_trajectories', all_trajectories);
save('blattner_exp_fit3/final_results.mat', 'final_results');

% Create a text file to save results
fileID = fopen('blattner_exp_fit3/detailed_results.txt', 'w');

% Write header
fprintf(fileID, 'Optimization Results Summary\n');
fprintf(fileID, '===========================\n\n');

% For each run
for i = 1:num_starts
    fprintf(fileID, 'Run %d:\n', i);
    fprintf(fileID, 'Parameters: s_bc=%.5f, s_cp=%.4f, s_A=%.5f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3));
    
    % Calculate individual errors for this solution
    [t_run, sol_run] = euler(@(t,y) model(t,y,solutions(i).params), [0, 24*100], [0,600,15.5], 0.01);
    
    % Extract last 36 hours and calculate errors (similar to your existing error calculation)
    csf_last_36hours = sol_run(233600:237200,2);
    
    time_indices = 1:200:3601;
    csf_model = csf_last_36hours(time_indices);
    
    % Ensure column vectors
    csf_model1 = csf_model(:);
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
fprintf(fileID, 'Parameters: s_bc=%.5f, s_cp=%.4f, s_A=%.5f\n', ...
    best_params(1), best_params(2), best_params(3));
fprintf(fileID, 'Objective value: %.6f\n', best_fval);

fclose(fileID);

% Save loss curves for each run
for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = arrayfun(@(idx) objective_function_fmincon(loss_curve(idx,:), csf_conc_exp1), iterations);
    loss_data = [iterations(:), loss_values(:)];
    writematrix(loss_data, sprintf('blattner_exp_fit3/loss_curve_run_%d.csv', i));
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
saveas(gcf, 'blattner_exp_fit3/loss_curves.png');

% Plot parameter space
plot_parameter_space(solutions, num_starts);

% Print final best results
fprintf('\nBest solution found:\n');
fprintf('Parameters: s_bc=%.5f, s_cp=%.5g, s_A=%.5f\n', ...
    best_params(1), best_params(2), best_params(3));
fprintf('Objective value: %.6f\n', best_fval);
