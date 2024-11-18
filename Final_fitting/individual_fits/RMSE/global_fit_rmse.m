global optimization_history;
optimization_history = [];

% Set paths to your plasma data files
csf_data_file1 = '/home/satyam/Documents/combined_fitting_data/final_fitting/blattner2020_csf42.csv';
csf_data_file2 = '/home/satyam/Documents/combined_fitting_data/final_fitting/lucey2018_csf42.csv';
csf_data_file3 = '/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_csf42.csv';
plasma_data_file1 = '/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_plasma42.csv';

% Read plasma data from both files
csf_data1 = readtable(csf_data_file1);
csf_data2 = readtable(csf_data_file2);
csf_data3 = readtable(csf_data_file3);
plasma_data1 = readtable(plasma_data_file1);

% Extract data from both plasma files
time_exp1 = csf_data1.Time;
time_exp2 = csf_data2.Time;
time_exp3 = csf_data3.Time;
csf_conc_exp1 = csf_data1.Conc;
csf_conc_exp2 = csf_data2.Conc;
csf_conc_exp3 = csf_data3.Conc;
plasma_conc_exp1 = plasma_data1.Conc;
csf_lsd1 = csf_data1.LSD;
csf_lsd2 = csf_data2.LSD;
csf_lsd3 = csf_data3.LSD;
csf_usd1 = csf_data1.USD;
csf_usd2 = csf_data2.USD;
csf_usd3 = csf_data3.USD;
plasma_lsd1 = plasma_data1.LSD;
plasma_usd1 = plasma_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;
csf_std2 = (csf_usd2 - csf_lsd2)/2;
csf_std3 = (csf_usd3 - csf_lsd3)/2;
plasma_std1 = (plasma_usd1 - plasma_lsd1)/2;

learning_rate = 0.01;

% Updated model function to optimize only sigma_bp (a), sigma_cp (b), and rbc (a12_wake)
function dydt_n = model(t, y, params)
    r_cp = params(1);
    sigma_bp = params(2);
    sigma_cp = params(3);
    sigma_p = params(4);

    % Derived parameters
    A = 13;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_p = 0.28;
    r_bp = (r_bc*(1-133*r_cp))/(133*r_cp);

    % Switch between sleep and wake states
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_cp * r_bc * (1 - sw_cycle)) * y(1);
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

% Modified objective function with error handling
function [total_error] = objective_function_fmincon(params, exp_csf1, exp_csf2, exp_csf3, exp_plasma1)
    try
        % Check for invalid parameter values
        if any(isnan(params)) || any(isinf(params)) || any(params < 0)
            total_error = 1e10;  % Return large value for invalid parameters
            return;
        end
        
        [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
        
        % Check for NaN or Inf in solution
        if any(isnan(sol(:))) || any(isinf(sol(:)))
            total_error = 1e10;
            return;
        end
        
        % Extract last 36 hours of data
        csf_last_36hours_data = sol(233600:237200,2);
        plasma_last_36hours_data = sol(233600:237200,3);
        
        % Select data points corresponding to experimental time points
        time_indices = 1:200:3601;
        csf_model = csf_last_36hours_data(time_indices);
        plasma_model = plasma_last_36hours_data(time_indices);
        
        % Calculate average of first 7 points
        avg_C_model = mean(csf_model(1:7));
        avg_P_model = mean(plasma_model(1:7));
        
        % Check for zero averages to prevent division by zero
        if avg_C_model == 0 || avg_P_model == 0
            total_error = 1e10;
            return;
        end
        
        % Normalize the model data
        normalised_C_model = (csf_model/avg_C_model)*100;
        normalised_P_model = (plasma_model/avg_P_model)*100;
        
        % Calculate error for each dataset
        error_C1 = sqrt(mean((normalised_C_model - exp_csf1).^2));
        error_C2 = sqrt(mean((normalised_C_model - exp_csf2).^2));
        error_C3 = sqrt(mean((normalised_C_model - exp_csf3).^2));
        error_P1 = sqrt(mean((normalised_P_model - exp_plasma1).^2));
        
        % Total error is the sum of errors from all datasets
        total_error = (error_C1 + error_C2 + error_C3)/3 + error_P1;
        
        % Final check for invalid result
        if isnan(total_error) || isinf(total_error)
            total_error = 1e10;
        end
    catch ME
        % If any error occurs, return a large value
        total_error = 1e10;
    end
end

% Modified sequential optimization function with better error handling
function [best_params, best_fval] = sequential_multistart_optimization(num_starts, initial_guess, lb, ub, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1)
    % Initialize storage for best results
    best_params = [];
    best_fval = Inf;
    
    % Setup optimization options with more robust settings
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'MaxIterations', 2000, ...
        'MaxFunctionEvaluations', 10000, ...
        'Algorithm', 'interior-point', ...
        'FiniteDifferenceStepSize', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6, ...
        'ScaleProblem', 'obj-and-constr', ...
        'TypicalX', initial_guess);
    
    % Create output directory
    if ~exist('global_optimization_results_rmse', 'dir')
        mkdir('global_optimization_results_rmse');
    end
    
    % Run optimizations sequentially with error handling
    for run = 1:num_starts
        fprintf('\nStarting optimization run %d/%d\n', run, num_starts);
        
        % Generate initial point
        if run == 1
            x0 = initial_guess;
        else
            % Random initial point with bounds checking
            max_attempts = 10;
            attempt = 1;
            valid_point_found = false;
            
            while attempt <= max_attempts && ~valid_point_found
                x0 = lb + rand(size(lb)) .* (ub - lb);
                % Test if initial point gives valid objective
                test_obj = objective_function_fmincon(x0, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1);
                if test_obj < 1e10
                    valid_point_found = true;
                end
                attempt = attempt + 1;
            end
            
            if ~valid_point_found
                fprintf('Warning: Could not find valid initial point for run %d\n', run);
                continue;
            end
        end
        
        try
            % Create objective function for this run
            this_run_objective = @(params) objective_function_wrapper(params, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1, run);
            
            % Run optimization with error catching
            [params, fval] = fmincon(this_run_objective, x0, [], [], [], [], lb, ub, [], options);
            
            % Check if optimization was successful
            if ~isempty(params) && ~isnan(fval) && ~isinf(fval) && fval < best_fval
                best_fval = fval;
                best_params = params;
                fprintf('New best solution found in run %d with loss: %f\n', run, fval);
            end
            
            % Save run results
            save(sprintf('global_optimization_results_rmse/run_%d_results.mat', run), 'params', 'fval', 'x0');
            
        catch ME
            fprintf('Error in optimization run %d: %s\n', run, ME.message);
            continue;
        end
    end
    
    % Check if any valid solution was found
    if isempty(best_params)
        error('No valid solution found in any optimization run');
    end
    
    % Save final results
    save('global_optimization_results_rmse/final_results.mat', 'best_params', 'best_fval');
end

% Modified objective function wrapper for sequential runs
function [total_error] = objective_function_wrapper(params, exp_csf1, exp_csf2, exp_csf3, exp_plasma1, run_number)
    % Calculate actual error
    total_error = objective_function_fmincon(params, exp_csf1, exp_csf2, exp_csf3, exp_plasma1);
    
    % Create persistent variables to track state
    persistent iteration_count;
    
    % Initialize on first call of each run
    if isempty(iteration_count)
        iteration_count = 1;
        
        % Create new file with headers
        filename = sprintf('global_optimization_results_rmse/run_%d_history.txt', run_number);
        fid = fopen(filename, 'w');
        fprintf(fid, 'Iteration\tLoss\tr_cp\tsigma_bp\tsigma_cp\tsigma_p\n');
        fprintf(fid, '%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n', ...
            iteration_count, total_error, params(1), params(2), params(3), params(4));
        fclose(fid);
    else
        % Increment iteration and log
        iteration_count = iteration_count + 1;
        
        filename = sprintf('global_optimization_results_rmse/run_%d_history.txt', run_number);
        fid = fopen(filename, 'a');
        fprintf(fid, '%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n', ...
            iteration_count, total_error, params(1), params(2), params(3), params(4));
        fclose(fid);
    end
    
    % Reset iteration count at the end of each run
    if total_error < 1e-6 || iteration_count >= 2000
        iteration_count = [];
    end
end

% Usage example:
num_starts = 20;
initial_guess = [0.005, 2, 2, 2];
lb = [0, 1, 1, 1];
ub = [0.0074, 10, 10, 10];

[best_params, best_fval] = sequential_multistart_optimization(num_starts, ...
    initial_guess, lb, ub, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1);

% Run the model for 100 days with optimized parameters
[t_100days, sol_100days] = euler(@(t,y) model(t,y,best_params), [0, 24*100], [0,600,15.5], 0.01);

% Save the simulation results using high precision
fileID = fopen('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values.csv', 'w');
fprintf(fileID, 'Time,Brain,CSF,Plasma\n');
for i = 1:length(t_100days)
    fprintf(fileID, '%.2f,%.15f,%.15f,%.15f\n', t_100days(i), sol_100days(i, 1), sol_100days(i, 2), sol_100days(i,3));
end
fclose(fileID);

% Print optimized parameters
fprintf('\nOptimized Parameters for Three Datasets:\n');
fprintf('r_cp: %f\n', best_params(1));
fprintf('sigma_bp: %f\n', best_params(2));
fprintf('sigma_cp: %f\n', best_params(3));
fprintf('sigma_p: %f\n', best_params(4));

% Run model with optimized parameters
[t, y] = euler(@(t, y) model(t, y, best_params), [0, 24*100], [0, 600, 15.5], 0.01);

% Load saved results for comparison
saved_data = readmatrix('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values.csv');
t_saved = saved_data(:, 1);
y_saved = saved_data(:, 2:end);

% Plot comparison
figure;
subplot(3, 1, 1);
hold on;
plot(t, y(:, 1), 'r', 'LineWidth', 2.0, 'DisplayName', 'Brain');
legend('show');
xlabel('Time (hr)');
ylabel('Brain Concentration');
title('Comparison of Model Simulations');
xlim([2336, 2386]);
xticks(2336:2:2386); % Set x-tick positions
xticklabels(0:2:50); % Set x-tick labels to start from 0
xline(2352, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2360, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2376, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2384, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;

subplot(3, 1, 2);
hold on;
plot(t, y(:, 2), 'r', 'LineWidth', 2.0, 'DisplayName', 'CSF');
legend('show');
xlabel('Time (hr)');
ylabel('CSF Concentration');
xlim([2336, 2386]);
xticks(2336:2:2386); % Set x-tick positions
xticklabels(0:2:50); % Set x-tick labels to start from 0
xline(2352, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2360, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2376, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2384, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;

subplot(3, 1, 3);
hold on;
plot(t, y(:, 3), 'r', 'LineWidth', 2.0, 'DisplayName', 'Plasma');
legend('show');
xlabel('Time (hr)');
ylabel('Plasma Concentration');
xlim([2336, 2386]);
xticks(2336:2:2386); % Set x-tick positions
xticklabels(0:2:50); % Set x-tick labels to start from 0
xline(2352, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2360, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2376, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2384, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;

cdata_12hours = sol_100days(233600:234800, 2);
cdata_36hours = sol_100days(233600:237200, 2);
pdata_12hours = sol_100days(233600:234800, 3);
pdata_36hours = sol_100days(233600:237200, 3);
csf_mean = mean(cdata_12hours);
plasma_mean = mean(pdata_12hours);
normalise_C = (cdata_36hours/ csf_mean)*100;
normalise_P = (pdata_36hours / plasma_mean)*100;
disp(length(normalise_C));
disp(length(normalise_P));

cdata_12hours_saved = y_saved(233600:234800, 2);
cdata_36hours_saved = y_saved(233600:237200, 2);
pdata_12hours_saved = y_saved(233600:234800, 3);
pdata_36hours_saved = y_saved(233600:237200, 3);
csf_mean_saved = mean(cdata_12hours_saved);
plasma_mean_saved = mean(pdata_12hours_saved);
normalise_C_saved = (cdata_36hours_saved / csf_mean_saved)*100;
normalise_P_saved = (pdata_36hours_saved / plasma_mean_saved)*100;
disp(length(normalise_C_saved));
disp(length(normalise_P_saved));

experimental_csf1 = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/blattner2020_csf42.csv');
experimental_csf2 = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/lucey2018_csf42.csv');
experimental_csf3 = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_csf42.csv');
experimental_plasma1 = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_plasma42.csv');
experimental_csf1.Time = experimental_csf1.Time - 100;
experimental_csf2.Time = experimental_csf2.Time - 100;
experimental_csf3.Time = experimental_csf3.Time - 100;
experimental_plasma1.Time = experimental_plasma1.Time - 100;
time_exp = experimental_plasma1.Time;
csf1_conc_exp = experimental_csf1.Conc;
csf2_conc_exp = experimental_csf2.Conc;
csf3_conc_exp = experimental_csf3.Conc;
plasma1_conc_exp = experimental_plasma1.Conc;

x = 1:3601;
ticks = 101:100:3601;
x = 1:3601;
disp(length(x));

% Define new ticks and labels
new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);

% Ensure ticks + 1 does not exceed the length of normalise_C
valid_indices = ticks + 1 <= length(normalise_C);
valid_ticks = ticks(valid_indices);

% Plot comparison CSF
figure;
subplot(3, 1, 1);
hold on;
plot(x, normalise_C, 'r', 'LineWidth', 2.0, 'HandleVisibility','off');
plot(x, normalise_C_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp/10, csf1_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
scatter(valid_ticks, normalise_C(valid_ticks + 1), 'r', 'DisplayName','Model Data');
axhline = refline(0, 0); 
axhline.LineStyle = '--'; 
axhline.Color = 'black';
axhline.HandleVisibility = 'off';
legend('show');
xlabel('Time (hr)');
ylabel('CSF Concentration');
title('Comparison of Model Simulations');
xlim([0, 3500]);

% Set x-ticks and x-tick labels for the first subplot
xticks(0:100:3500); % Adjust the interval as needed
xticklabels(0:1:36); % Renumber x-ticks from 0 to 36

hold off;

valid_indices = ticks + 1 <= length(normalise_C); % Ensure valid indices
valid_ticks = ticks(valid_indices);
subplot(3, 1, 2);
hold on;
plot(x, normalise_C, 'r', 'LineWidth', 2.0, 'HandleVisibility', 'off');
plot(x, normalise_C_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp/10, csf2_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Lucey2022 - CSF');
scatter(valid_ticks, normalise_C(valid_ticks + 1), 'r', 'DisplayName', 'Model Data');
axhline = refline(0, 0); 
axhline.LineStyle = '--'; 
axhline.Color = 'black';
axhline.HandleVisibility = 'off';
legend('show');
xlabel('Time (hr)');
ylabel('CSF Concentration');
xlim([0, 3500]);

% Set x-ticks and x-tick labels for the second subplot
xticks(0:100:3500); % Adjust the interval as needed
xticklabels(0:1:36); % Renumber x-ticks from 0 to 36

hold off;

% Plot comparison CSF
subplot(3, 1, 3);
hold on;
plot(x, normalise_C, 'r', 'LineWidth', 2.0, 'HandleVisibility','off');
plot(x, normalise_C_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp/10, csf3_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Liu2018 - CSF');
scatter(valid_ticks, normalise_C(valid_ticks + 1), 'r', 'DisplayName','Model Data');
axhline = refline(0, 0); 
axhline.LineStyle = '--'; 
axhline.Color = 'black';
axhline.HandleVisibility = 'off';
legend('show');
xlabel('Time (hr)');
ylabel('CSF Concentration');
title('Comparison of Model Simulations');
xlim([0, 3500]);

% Set x-ticks and x-tick labels for the first subplot
xticks(0:100:3500); % Adjust the interval as needed
xticklabels(0:1:36); % Renumber x-ticks from 0 to 36

hold off;

valid_indices = ticks + 1 <= length(normalise_P); % Ensure valid indices
valid_ticks = ticks(valid_indices);

figure;
hold on;
plot(x, normalise_P, 'r', 'LineWidth', 2.0, 'HandleVisibility', 'off');
plot(x, normalise_P_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp/10, plasma1_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Liu Data - Plasma');
scatter(valid_ticks, normalise_P(valid_ticks + 1), 'r', 'DisplayName', 'Model Data');
axhline = refline(0, 0); 
axhline.LineStyle = '--'; 
axhline.Color = 'black';
axhline.HandleVisibility = 'off';
legend('show');
xlabel('Time (hr)');
ylabel('Plasma Concentration');
xlim([0, 3500]);

% Set x-ticks and x-tick labels for the second subplot
xticks(0:100:3500); % Adjust the interval as needed
xticklabels(0:1:36); % Renumber x-ticks from 0 to 36
hold off;