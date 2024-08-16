% Clear the loss_values array at the beginning of the script
clear all;
clear global loss_values;
global loss_values;
loss_values = [];

% Read experimental data from CSV files
data_exp1 = readtable('/home/satyam/Documents/MATLAB/Model_fitting_SD/liu2022_csf.csv');
data_exp2 = readtable('/home/satyam/Documents/MATLAB/Model_fitting_SD/liu2022_plasma.csv');
time_exp = data_exp1.Time;
exp_csf = data_exp1.Conc;
exp_plasma = data_exp2.Conc;
csf_lsd = data_exp1.LSD;
csf_usd = data_exp1.USD;
plasma_lsd = data_exp2.LSD;
plasma_usd = data_exp2.USD;
csf_std = (csf_usd - csf_lsd)/2;
plasma_std = (plasma_usd - plasma_lsd)/2;

% Define the model function
function dydt_n = model(t, y, params)
    a = params(1);
    b = params(2);
    c = params(3);
    A_wake = params(4);
    A_sleep = 0.8*A_wake;
    a12_wake = params(5);
    a12_sleep = 2.5*a12_wake;
    a13_wake = params(6);
    a13_sleep = a*a13_wake;
    a23_sleep = 0.06601;
    a23_wake = a23_sleep/b;
    k_wake = 0.231049;
    k_sleep = c*k_wake;

    sw_cycle = (mod(t,24) >= 8 && mod(t,24) < 24);
    dydt_n = zeros(3, 1);
    dydt_n(1) = A_wake*sw_cycle + (A_sleep*(1 - sw_cycle)) - (a13_wake*sw_cycle + a13_sleep*(1-sw_cycle) + a12_wake*sw_cycle + a12_sleep*(1 - sw_cycle))*y(1);
    dydt_n(2) = (a12_wake*sw_cycle + a12_sleep*(1-sw_cycle))*y(1) - (a23_wake * sw_cycle + a23_sleep*(1-sw_cycle))*y(2);
    dydt_n(3) = (a23_wake * sw_cycle + a23_sleep*(1-sw_cycle))*y(2) + (a13_wake*sw_cycle + a13_sleep*(1-sw_cycle))*y(1) - (k_wake*sw_cycle + k_sleep*(1 - sw_cycle))*y(3);
end

% Define the Euler-Maruyama method
function [t, w] = euler(F, endpoints, parameters, ts)
    if length(endpoints) == 2
        h = ts; % deltat (seconds)
        if mod((endpoints(2) - endpoints(1))/h,1)
            error('Step size does not divide into endpoints evenly');
        end
        t = endpoints(1):h:endpoints(2);
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints;
    end

    w = zeros(length(t), length(parameters));
    w(1,:) = parameters;

    for k = 1:length(t)-1
        w(k+1,:) = w(k,:) + feval(F,t(k),w(k,:))'*h;
    end

    t = t(:);
end

% Define the objective function
function [total_error, is_valid] = objective_function(params, exp_csf, exp_plasma, csf_std, plasma_std)
    global loss_values;

    % Extract parameters
    a = params(1);
    b = params(2);
    c = params(3);
    A_wake = params(4);
    a12_wake = params(5);
    a13_wake = params(6);

    % Check constraints
    wake_constraint = 0 <= a12_wake / (a12_wake + a13_wake) && a12_wake / (a12_wake + a13_wake) <= 4.26 / b;
    sleep_constraint = 0 <= (2.5 * a12_wake) / (2.5 * a12_wake + a * a13_wake) && (2.5 * a12_wake) / (2.5 * a12_wake + a * a13_wake) <= 5.336;

    is_valid = wake_constraint && sleep_constraint;

    if ~is_valid
        total_error = Inf;
        return;
    end

    % Use Euler-Maruyama method instead of ode45
    [~, sol] = euler(@(t, y) model(t, y, params), [0, 24*100], [0, 600, 15.5], 0.01);
    
    csf_last_36hours_data = sol(236100:239700, 2);
    plasma_last_36hours_data = sol(236100:239700, 3);

    time_indices = 1:200:3601;
    csf_last_36hours = csf_last_36hours_data(time_indices);
    plasma_last_36hours = plasma_last_36hours_data(time_indices);

    avg_C = mean(csf_last_36hours(1:6));
    avg_B = mean(plasma_last_36hours(1:6));

    normalised_C = (csf_last_36hours / avg_C) * 100;
    normalised_B = (plasma_last_36hours / avg_B) * 100;

    csf_weights = 1./(csf_std.^2);
    plasma_weights = 1./(plasma_std.^2);

    error_C = sqrt(mean(csf_weights .* (normalised_C - exp_csf).^2));
    error_B = sqrt(mean(plasma_weights .* (normalised_B - exp_plasma).^2));

    total_error = error_C + error_B;
    
    % Update the global loss_values for this iteration
    loss_values(end+1) = total_error;
end

% Define initial guess for parameters
initial_guess = [3, 4, 4, 21, 0.3, 0.02];

% Define bounds for parameters
lb = [1, 1, 1, 15, 0.1, 0.01];
ub = [10, 10, 10, 24, 1, 0.1];

% Define the objective function for fmincon
function [total_error] = fmincon_objective(params, exp_csf, exp_plasma, csf_std, plasma_std)
    [total_error, is_valid] = objective_function(params, exp_csf, exp_plasma, csf_std, plasma_std);
    if ~is_valid
        total_error = Inf;
    end
end

% Define constraint function
function [c, ceq] = constraint_fun(params)
    a = params(1);
    b = params(2);
    a12_wake = params(5);
    a13_wake = params(6);
    
    c1 = a12_wake / (a12_wake + a13_wake) - 4.26 / b;
    c2 = -(2.5 * a12_wake) / (2.5 * a12_wake + a * a13_wake);
    c3 = (2.5 * a12_wake) / (2.5 * a12_wake + a * a13_wake) - 5.336;
    
    c = [c1; c2; c3];
    ceq = [];
end

% Set up optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 100, 'MaxFunctionEvaluations', 10000);

% Run optimization
[optimized_params, fval] = fmincon(@(params) fmincon_objective(params, exp_csf, exp_plasma, csf_std, plasma_std), ...
                                   initial_guess, [], [], [], [], lb, ub, @constraint_fun, options);

% Save optimized parameters using high precision
fileID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/optimized_params_model_v51_clearance.csv', 'w');
fprintf(fileID, '%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n', optimized_params);
fclose(fileID);

% Save the loss values
lossID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v51_loss_clearance.csv', 'w');
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
[t, sol] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 600, 15.5], 0.01);

% Save the simulation results using high precision
fileID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v51_output_clearance.csv', 'w');
fprintf(fileID, 'Time,Brain,CSF,Plasma\n');
for i = 1:length(t)
    fprintf(fileID, '%.2f,%.15f,%.15f,%.15f\n', t(i), sol(i, 1), sol(i, 2), sol(i, 3));
end
fclose(fileID);

% Print optimized parameters
fprintf('Optimized Parameters:\n');
fprintf('a: %f\n', optimized_params(1));
fprintf('b: %f\n', optimized_params(2));
fprintf('c: %f\n', optimized_params(3));
fprintf('A_wake: %f\n', optimized_params(4));
fprintf('a12_wake: %f\n', optimized_params(5));
fprintf('a13_wake: %f\n', optimized_params(6));

% Load optimized parameters
optimized_params = readmatrix('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/optimized_params_model_v51_clearance.csv');

% Run model with optimized parameters
[t, y] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 600, 15.5], 0.01);

% Load saved results for comparison
saved_data = readmatrix('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v51_output_clearance.csv');
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
normalise_C = (csf_data_36hours / csf_mean) * 100;
normalise_B = (plasma_data_36hours / plasma_mean) * 100;
disp(length(normalise_C));

csf_saved_12hours = y_saved(236100:237300, 2);
plasma_saved_12hours = y_saved(236100:237300, 3);
csf_saved_36hours = y_saved(236100:239700, 2);
plasma_saved_36hours = y_saved(236100:239700, 3);
csf_saved_mean = mean(csf_saved_12hours);
plasma_saved_mean = mean(plasma_saved_12hours);
normalise_C_saved = (csf_saved_36hours / csf_saved_mean) * 100;
normalise_B_saved = (plasma_saved_36hours / plasma_saved_mean) * 100;
disp(length(normalise_C_saved));

experimental_csf = readtable('/home/satyam/Documents/MATLAB/Model_fitting_SD/liu2022_csf_1000.csv');
experimental_plasma = readtable('/home/satyam/Documents/MATLAB/Model_fitting_SD/liu2022_plasma_1000.csv');
experimental_csf.Time = experimental_csf.Time - 100;
experimental_plasma.Time = experimental_plasma.Time - 100;
time_exp = experimental_csf.Time;
csf_conc_exp = experimental_csf.Conc;
plasma_conc_exp = experimental_plasma.Conc;

% Sample the normalized data at 2-hour intervals
time_indices = 1:200:3601; % This will give 19 points over 36 hours
normalise_C_sampled = normalise_C(time_indices);
normalise_B_sampled = normalise_B(time_indices);

x = 1:3601;
ticks = 101:100:3601;
x = 1:3601;
disp(length(x));

time_sampled = 0:200:3601;

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
errorbar(time_exp, csf_conc_exp, csf_conc_exp - csf_lsd, csf_usd - csf_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Experimental Data - CSF');
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
errorbar(time_exp, plasma_conc_exp, plasma_conc_exp - plasma_lsd, plasma_usd - plasma_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Experimental Data - Plasma');
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
