% Set paths to your plasma data files
csf_data_file1 = '/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_csf42.csv';
plasma_data_file1 = '/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_plasma42.csv';

% Read plasma data from both files
csf_data1 = readtable(csf_data_file1);
plasma_data1 = readtable(plasma_data_file1);

% Extract data from both plasma files
time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Conc;
plasma_conc_exp1 = plasma_data1.Conc;

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
function [t, w] = euler(F, endpoints, parameters, ts)
    if length(endpoints) == 2
        h = ts; %delta_t (seconds)
        total_time = endpoints(2) - endpoints(1);
        num_steps = floor(total_time / h);
        t = linspace(endpoints(1), endpoints(2), num_steps + 1); %Creating time vector
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints;
    end
    w = zeros(length(t), length(parameters));
    w(1,:) = parameters;
    for k = 1:length(t) - 1
        w(k+1,:) = w(k,:) + feval(F, t(k), w(k,:))' * h;
    end
    t = t(:);
end

function total_error = objective_function_fmincon(params, exp_csf1, exp_plasma1)
    try
        % Check if parameters are valid
        if any(isnan(params)) || any(isinf(params))
            total_error = Inf;
            return;
        end
        
        % Run simulation
        [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
        
        % Extract last 36 hours of data
        if size(sol, 1) >= 237200
            csf_last_36hours_data = sol(233600:237200,2);
            plasma_last_36hours_data = sol(233600:237200,3);
            
            % Select data points corresponding to experimental time points
            time_indices = 1:200:3601;
            if length(csf_last_36hours_data) >= max(time_indices)
                csf_model = csf_last_36hours_data(time_indices);
                plasma_model = plasma_last_36hours_data(time_indices);
                
                % Calculate average of first 7 points
                avg_C_model = mean(csf_model(1:7));
                avg_P_model = mean(plasma_model(1:7));

                % Normalize the model data
                if avg_C_model > 0
                    normalised_C_model = (csf_model/avg_C_model)*100;
                    normalised_P_model = (plasma_model/avg_P_model)*100;
                    
                    % Calculate error
                    if length(normalised_C_model) == length(exp_csf1)
                        C_error = sqrt(mean((normalised_C_model - exp_csf1).^2));
                        P_error = sqrt(mean((normalised_P_model - exp_plasma1).^2));
                        total_error = C_error + P_error;
                        return;
                    end
                end
            end
        end
        
        % If any conditions aren't met, return Inf
        total_error = Inf;
        
    catch
        total_error = Inf;
    end
end

% Main optimization function
function [best_params, best_fval] = sequential_multistart_optimization(num_starts, initial_guess, lb, ub, csf_conc_exp1, plasma_conc_exp1)
    % Initialize storage
    best_params = [];
    best_fval = Inf;
    
    % Setup optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'MaxIterations', 2000, ...
        'MaxFunctionEvaluations', 10000, ...
        'Algorithm', 'interior-point', ...
        'FiniteDifferenceStepSize', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-6);
    
    % Create output directory
    if ~exist('liu_optimization_results_rmse', 'dir')
        mkdir('liu_optimization_results_rmse');
    end
    
    % Run optimizations sequentially
    for run = 1:num_starts
        fprintf('\nStarting optimization run %d/%d\n', run, num_starts);
        
        % Generate initial point
        if run == 1
            x0 = initial_guess;
        else
            x0 = lb + rand(size(lb)) .* (ub - lb);
        end
        
        % Verify initial point is valid
        initial_error = objective_function_fmincon(x0, csf_conc_exp1, plasma_conc_exp1);
        if isinf(initial_error)
            fprintf('Skipping run %d due to invalid initial point\n', run);
            continue;
        end
        
        % Create objective function
        this_run_objective = @(params) objective_function_wrapper(params, csf_conc_exp1, plasma_conc_exp1, run);
        
        try
            % Run optimization
            [params, fval] = fmincon(this_run_objective, x0, [], [], [], [], lb, ub, [], options);
            
            % Update best result
            if fval < best_fval
                best_fval = fval;
                best_params = params;
                fprintf('New best solution found in run %d with loss: %f\n', run, fval);
            end
            
            % Save run results
            save(sprintf('liu_optimization_results_rmse/run_%d_results.mat', run), 'params', 'fval', 'x0');
            
        catch ME
            fprintf('Error in run %d: %s\n', run, ME.message);
            continue;
        end
    end
    
    % Save final results
    if ~isempty(best_params)
        save('liu_optimization_results_rmse/final_results.mat', 'best_params', 'best_fval');
    else
        error('No valid solution found in any optimization run');
    end
end

% Objective function wrapper
function total_error = objective_function_wrapper(params, exp_csf1, exp_plasma1, run_number)
    % Calculate error
    total_error = objective_function_fmincon(params, exp_csf1, exp_plasma1);
    
    % Create persistent counter
    persistent iteration_count;
    
    % Initialize or increment counter
    if isempty(iteration_count)
        iteration_count = 1;
        filename = sprintf('liu_optimization_results_rmse/run_%d_history.txt', run_number);
        fid = fopen(filename, 'w');
        fprintf(fid, 'Iteration\tLoss\tr_cp\tsigma_bp\tsigma_cp\tsigma_p\n');
    else
        iteration_count = iteration_count + 1;
        filename = sprintf('liu_optimization_results_rmse/run_%d_history.txt', run_number);
        fid = fopen(filename, 'a');
    end
    
    % Log results
    fprintf(fid, '%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n', ...
        iteration_count, total_error, params(1), params(2), params(3), params(4));
    fclose(fid);
    
    % Reset counter if needed
    if total_error < 1e-6 || iteration_count >= 2000
        iteration_count = [];
    end
end

% Usage example:
num_starts = 20;
initial_guess = [0.005, 3, 5, 4];
lb = [0, 1, 1, 1];
ub = [0.0074, 10, 10, 10];

[best_params, best_fval] = sequential_multistart_optimization(num_starts, initial_guess, lb, ub, csf_conc_exp1, plasma_conc_exp1);

% Run the model for 100 days with optimized parameters
[t_100days, sol_100days] = euler(@(t,y) model(t,y,best_params), [0, 24*100], [0,600,15.5], 0.01);

% Save the simulation results using high precision
fileID = fopen('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values_liu2022.csv', 'w');
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
[t, y] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 600, 15.5], 0.01);

% Load saved results for comparison
saved_data = readmatrix('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values_liu2022.csv');
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
normalise_P = (pdata_36hours/ plasma_mean)*100;
disp(length(normalise_C));

cdata_12hours_saved = y_saved(233600:234800, 2);
cdata_36hours_saved = y_saved(233600:237200, 2);
pdata_12hours_saved = y_saved(233600:234800, 3);
pdata_36hours_saved = y_saved(233600:237200, 3);
csf_mean_saved = mean(cdata_12hours_saved);
plasma_mean_saved = mean(pdata_12hours_saved);
normalise_C_saved = (cdata_36hours_saved / csf_mean_saved)*100;
normalise_P_saved = (pdata_36hours_saved / plasma_mean_saved)*100;
disp(length(normalise_C_saved));

experimental_csf1 = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_csf42.csv');
experimental_plasma1 = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_plasma42.csv');
experimental_csf1.Time = experimental_csf1.Time - 100;
experimental_plasma1.Time = experimental_plasma1.Time - 100;
time_exp = experimental_csf1.Time;
csf1_conc_exp = experimental_csf1.Conc;
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

figure;
hold on;
plot(x, normalise_P, 'r', 'LineWidth', 2.0, 'HandleVisibility','off');
plot(x, normalise_P_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp/10, plasma1_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
scatter(valid_ticks, normalise_P(valid_ticks + 1), 'r', 'DisplayName','Model Data');
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