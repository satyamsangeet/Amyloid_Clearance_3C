% Clear the loss_values array at the beginning of the script
clear all;
clear global loss_values;
global loss_values;
loss_values = [];

% Read experimental data from CSV files
data_exp1 = csvread('/home/satyam/Documents/huang_csf42.csv', 1, 0);
data_exp2 = csvread('/home/satyam/Documents/huang_plasma42.csv', 1, 0);
time_exp = data_exp1(:, 1);
exp_csf = data_exp1(:, 2);
exp_plasma = data_exp2(:, 2);

% Define the model function
function dydt_n = model(t, y, params)
    a12_wake = params(1);
    a13_wake = params(2);
    a23_wake = params(3);
    A_wake = params(4);
    A_sleep = params(5);
    a = params(6);
    k = params(7);
    
    sw_cycle = (mod(t,24) >= 8 && mod(t,24) < 24);
    dydt_n = zeros(3, 1);
    dydt_n(1) = A_wake*sw_cycle + (A_sleep*(1 - sw_cycle)) - (a13_wake*sw_cycle + a*a13_wake*(1-sw_cycle) + a12_wake*sw_cycle + a*a12_wake*(1 - sw_cycle))*y(1);
    dydt_n(2) = (a12_wake*sw_cycle + a*a12_wake*(1-sw_cycle))*y(1) - (a23_wake*sw_cycle + a*a23_wake*(1-sw_cycle))*y(2);
    dydt_n(3) = (a23_wake*sw_cycle + a*a23_wake*(1-sw_cycle))*y(2) + (a13_wake*sw_cycle + a*a13_wake*(1-sw_cycle))*y(1) - k*y(3);
end

% Define the Euler-Maruyama method
function [t, w] = euler(F, endpoints, parameters, ts)
    if length(endpoints) == 2
        h = ts; % deltat (seconds)
        total_time = endpoints(2) - endpoints(1);
        num_steps = floor(total_time / h); % Calculate the number of steps
        t = linspace(endpoints(1), endpoints(2), num_steps + 1); % Create a time vector
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints;
    end

    w = zeros(length(t), length(parameters)); % initialize output array
    w(1,:) = parameters;

    for k = 1:length(t) - 1
        w(k+1,:) = w(k,:) + feval(F, t(k), w(k,:))' * h; % Euler-Maruyama method
    end

    t = t(:);
end

% Define the objective function
function total_error = objective_function(params, exp_csf, exp_plasma)
    global loss_values;

    % Use Euler-Maruyama method instead of ode45
    [~, sol] = euler(@(t, y) model(t, y, params), [0, 24*100], [0, 0, 0], 0.01);
    
    csf_last_36hours_data = sol(236100:239700, 2);
    plasma_last_36hours_data = sol(236100:239700, 3);

    time_indices = 1:36;
    csf_last_36hours = csf_last_36hours_data(time_indices);
    plasma_last_36hours = plasma_last_36hours_data(time_indices);

    avg_C = mean(csf_last_36hours(1:12));
    avg_B = mean(plasma_last_36hours(1:12));

    normalised_C = (csf_last_36hours - avg_C) / avg_C;
    normalised_B = (plasma_last_36hours - avg_B) / avg_B;

    error_C = sqrt(mean((normalised_C - exp_csf).^2));
    error_B = sqrt(mean((normalised_B - exp_plasma).^2));
    total_error = error_C + error_B;
    
    % Update the global loss_values for this iteration
    loss_values(end+1) = total_error;
end

% Define initial guess for parameters
initial_guess = [0.5, 0.72, 0.3, 10, 7.2, 1.01, 0.25];

% Define bounds for parameters
lb = [0.3, 0.6, 0.2, 10, 7, 1.01, 0.2];
ub = [0.6, 0.9, 0.5, 11, 9, 1.05, 0.5];

% Custom optimization loop using Euler-Maruyama method
params = initial_guess;
learning_rate = 0.01;
num_iterations = 2000; % Set the correct number of iterations

for iter = 1:num_iterations
    current_error = objective_function(params, exp_csf, exp_plasma);
    
    % Calculate gradient using finite differences
    grad = zeros(size(params));
    for i = 1:length(params)
        params_temp = params;
        params_temp(i) = params_temp(i) + learning_rate;
        grad(i) = (objective_function(params_temp, exp_csf, exp_plasma) - current_error) / learning_rate;
    end
    
    % Update parameters
    params = params - learning_rate * grad;
    
    % Ensure parameters are within bounds
    params = max(min(params, ub), lb);
    
    fprintf('Iteration %d, Loss: %f\n', iter, current_error);
end

% Get optimized parameters
optimized_params = params;

% Save optimized parameters using high precision
fileID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/optimized_params.csv', 'w');
fprintf(fileID, '%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n', optimized_params);
fclose(fileID);

% Save the loss values
lossID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v43_loss.csv', 'w');
if lossID == -1
    disp('Error: Unable to open file for writing.');
else
    fprintf(lossID, 'Iteration,Loss\n');
    for i = 1:length(loss_values)
        try
            fprintf(lossID, '%d,%f\n', i, loss_values(i));
            disp(['Loss value at iteration ', num2str(i), ' saved successfully.']);
        catch
            fprintf('Error occurred while saving loss value at iteration %d\n', i);
        end
    end
    fclose(lossID);
    disp('Loss values saved successfully.');
end

% Solve the model with optimized parameters and save the results
[t, sol] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 0, 0], 0.01);

% Save the simulation results using high precision
fileID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v43_output.csv', 'w');
fprintf(fileID, 'Time,Brain,CSF,Plasma\n');
for i = 1:length(t)
    fprintf(fileID, '%.2f,%.15f,%.15f,%.15f\n', t(i), sol(i, 1), sol(i, 2), sol(i, 3));
end
fclose(fileID);

% Print optimized parameters
fprintf('Optimized Parameters:\n');
fprintf('A_wake: %f\n', optimized_params(4));
fprintf('A_sleep: %f\n', optimized_params(5));
fprintf('a12_wake: %f\n', optimized_params(1));
fprintf('k: %f\n', optimized_params(7));
fprintf('a: %f\n', optimized_params(6));
fprintf('a13_wake: %f\n', optimized_params(2));
fprintf('a23_wake: %f\n', optimized_params(3));

% Load optimized parameters
optimized_params = readmatrix('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/optimized_params.csv');

% Run model with optimized parameters
[t, y] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 0, 0], 0.01);

% Load saved results for comparison
saved_data = readmatrix('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v43_output.csv');
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


% Normalize and plot CSF and Plasma data
csf_data_12hours = y(236100:237300, 2);
plasma_data_12hours = y(236100:237300, 3);
csf_data_36hours = y(236100:239700, 2);
plasma_data_36hours = y(236100:239700, 3);
csf_mean = mean(csf_data_12hours);
plasma_mean = mean(plasma_data_12hours);
normalise_C = (csf_data_36hours - csf_mean) / csf_mean;
normalise_B = (plasma_data_36hours - plasma_mean) / plasma_mean;
disp(length(normalise_C));

csf_saved_12hours = y_saved(236100:237300, 2);
plasma_saved_12hours = y_saved(236100:237300, 3);
csf_saved_36hours = y_saved(236100:239700, 2);
plasma_saved_36hours = y_saved(236100:239700, 3);
csf_saved_mean = mean(csf_saved_12hours);
plasma_saved_mean = mean(plasma_saved_12hours);
normalise_C_saved = (csf_saved_36hours - csf_saved_mean) / csf_saved_mean;
normalise_B_saved = (plasma_saved_36hours - plasma_saved_mean) / plasma_saved_mean;
disp(length(normalise_C_saved));

experimental_csf = readtable('/home/satyam/Documents/huang_csf42_1000.csv');
experimental_plasma = readtable('/home/satyam/Documents/huang_plasma42_1000.csv');
experimental_csf.Time = experimental_csf.Time - 100;
experimental_plasma.Time = experimental_plasma.Time - 100;
time_exp = experimental_csf.Time;
csf_conc_exp = experimental_csf.Conc;
plasma_conc_exp = experimental_plasma.Conc;

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

% Plot comparison
figure;
subplot(2, 1, 1);
hold on;
plot(x, normalise_C, 'r', 'LineWidth', 2.0, 'HandleVisibility','off');
plot(x, normalise_C_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp, csf_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Experimental Data - CSF');
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

valid_indices = ticks + 1 <= length(normalise_B); % Ensure valid indices
valid_ticks = ticks(valid_indices);
subplot(2, 1, 2);
hold on;
plot(x, normalise_B, 'r', 'LineWidth', 2.0, 'HandleVisibility', 'off');
plot(x, normalise_B_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp, plasma_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Experimental Data - Plasma');
scatter(valid_ticks, normalise_B(valid_ticks + 1), 'r', 'DisplayName', 'Model Data');
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
