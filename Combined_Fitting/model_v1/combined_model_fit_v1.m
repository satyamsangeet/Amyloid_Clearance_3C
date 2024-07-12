% Set paths to your plasma data files
plasma_data_file1 = '/home/satyam/Documents/combined_fitting_data/huang_plasma42_new.csv';
plasma_data_file2 = '/home/satyam/Documents/combined_fitting_data/liu_plasma42.csv';

% Read plasma data from both files
plasma_data1 = readtable(plasma_data_file1);
plasma_data2 = readtable(plasma_data_file2);

% Extract data from both plasma files
time_exp1 = plasma_data1.Time;
plasma_conc_exp1 = plasma_data1.Conc;
time_exp2 = plasma_data2.Time;
plasma_conc_exp2 = plasma_data2.Conc;

% Initial parameters for model fitting
initial_guess = [0.5, 10, 7.2, 1.02, 0.25];
lb = [0.1, 1, 1, 1.01, 0.1];
ub = [1, 20, 10, 1.05, 1];
learning_rate = 0.01;
num_iterations = 2000;

% Define the model function
function dydt_n = model(t, y, params)
    a12_wake = params(1);
    A_wake = params(2);
    A_sleep = params(3);
    a = params(4);
    k = params(5);
    
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    dydt_n = zeros(2, 1);
    dydt_n(1) = A_wake * sw_cycle + A_sleep * (1 - sw_cycle) - (a12_wake * sw_cycle + a * a12_wake * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (a12_wake * sw_cycle + a * a12_wake * (1 - sw_cycle)) * y(1) - k * y(2);
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

% Defining Objective Function
function total_error = objective_function(params, exp_plasma1, exp_plasma2)
    global loss_values;
    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,0], 0.01);
    
    % Extract last 36 hours of data
    plasma_last_36hours_data = sol(236100:239700,2);
    
    % Select data points corresponding to experimental time points (every 2 hours)
    time_indices = 1:200:3601; % This will give 19 points over 36 hours (0 to 36, step 2)
    plasma_model = plasma_last_36hours_data(time_indices);
    
    % Calculate average of first 7 points for model data
    avg_P_model = mean(plasma_model(1:7));
    
    % Normalize the model data
    normalised_P_model = (plasma_model/avg_P_model)*100;
    
    % Calculate error for dataset 1
    avg_P_exp1 = mean(exp_plasma1(1:7));
    normalised_P_exp1 = (exp_plasma1/avg_P_exp1)*100;
    error_P1 = sqrt(mean((normalised_P_model - normalised_P_exp1).^2));
    
    % Calculate error for dataset 2
    avg_P_exp2 = mean(exp_plasma2(1:7));
    normalised_P_exp2 = (exp_plasma2/avg_P_exp2)*100;
    error_P2 = sqrt(mean((normalised_P_model - normalised_P_exp2).^2));
    
    % Total error is the sum of errors from both datasets
    total_error = error_P1 + error_P2;
    
    loss_values(end+1) = total_error;
end

% Main optimization loop for both datasets combined
params = initial_guess;
global loss_values;
loss_values = [];

for iter = 1:num_iterations
    current_error = objective_function(params, plasma_conc_exp1, plasma_conc_exp2);
    
    % Calculate gradient using finite differences
    grad = zeros(size(params));
    for i = 1:length(params)
        params_temp = params;
        params_temp(i) = params_temp(i) + learning_rate;
        grad(i) = (objective_function(params_temp, plasma_conc_exp1, plasma_conc_exp2) - current_error) / learning_rate;
    end
    
    % Update parameters
    params = params - learning_rate * grad;
    
    % Ensure parameters are within bounds
    params = max(min(params, ub), lb);
    
    fprintf('Iteration %d, Loss: %f\n', iter, current_error);
end

% Optimized parameters
optimized_params = params;

% Run the model for 100 days with optimized parameters
[t_100days, sol_100days] = euler(@(t,y) model(t,y,optimized_params), [0, 24*100], [0,0], 0.01);

% Save the simulation results using high precision
fileID = fopen('/home/satyam/Documents/combined_fitting_data/combined_model_v1_output.csv', 'w');
fprintf(fileID, 'Time,Brain,Plasma\n');
for i = 1:length(t)
    fprintf(fileID, '%.2f,%.15f,%.15f\n', t_100days(i), sol_100days(i, 1), sol_100days(i, 2));
end
fclose(fileID);

% Print optimized parameters
fprintf('\nOptimized Parameters for Both Datasets:\n');
fprintf('a12_wake: %f\n', optimized_params(1));
fprintf('A_wake: %f\n', optimized_params(2));
fprintf('A_sleep: %f\n', optimized_params(3));
fprintf('a: %f\n', optimized_params(4));
fprintf('k: %f\n', optimized_params(5));

% Run model with optimized parameters
[t, y] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 0], 0.01);

% Load saved results for comparison
saved_data = readmatrix('/home/satyam/Documents/combined_fitting_data/combined_model_v1_output.csv');
t_saved = saved_data(:, 1);
y_saved = saved_data(:, 2:end);

% Plot comparison
figure;
subplot(2, 1, 1);
hold on;
plot(t, y(:, 1), 'r', 'LineWidth', 2.0, 'DisplayName', 'CSF');
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

subplot(2, 1, 2);
hold on;
plot(t, y(:, 2), 'r', 'LineWidth', 2.0, 'DisplayName', 'Plasma');
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

pdata_12hours = sol_100days(236100:237300, 2);
pdata_36hours = sol_100days(236100:239700, 2);
plasma_mean = mean(pdata_12hours);
normalise_P = (pdata_36hours / plasma_mean)*100;
disp(length(normalise_P));

pdata_12hours_saved = y_saved(236100:237300, 2);
pdata_36hours_saved = y_saved(236100:239700, 2);
plasma_mean_saved = mean(pdata_12hours);
normalise_P_saved = (pdata_36hours / plasma_mean)*100;
disp(length(normalise_P_saved));

experimental_plasma1 = readtable('/home/satyam/Documents/combined_fitting_data/huang_plasma42_1000.csv');
experimental_plasma2 = readtable('/home/satyam/Documents/combined_fitting_data/liu_plasma42_1000.csv');
experimental_plasma1.Time = experimental_plasma1.Time - 100;
experimental_plasma2.Time = experimental_plasma2.Time - 100;
time_exp = experimental_plasma1.Time;
plasma1_conc_exp = experimental_plasma1.Conc;
plasma2_conc_exp = experimental_plasma2.Conc;

x = 1:3601;
ticks = 101:100:3601;
x = 1:3601;
disp(length(x));

% Define new ticks and labels
new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);

% Ensure ticks + 1 does not exceed the length of normalise_C
valid_indices = ticks + 1 <= length(normalise_P);
valid_ticks = ticks(valid_indices);

% Plot comparison
figure;
subplot(2, 1, 1);
hold on;
plot(x, normalise_P, 'r', 'LineWidth', 2.0, 'HandleVisibility','off');
plot(x, normalise_P_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp, plasma1_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Huang Data - Plasma');
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

valid_indices = ticks + 1 <= length(normalise_P); % Ensure valid indices
valid_ticks = ticks(valid_indices);
subplot(2, 1, 2);
hold on;
plot(x, normalise_P, 'r', 'LineWidth', 2.0, 'HandleVisibility', 'off');
plot(x, normalise_P_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp, plasma2_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Liu Data - Plasma');
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
