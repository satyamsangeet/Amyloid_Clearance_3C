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

%initial_guess = [1.1, 1.2, 1.3, 0.3, 0.02, 22];
%lb = [1, 1, 1, 0.1, 0.01, 15];
%ub = [3, 3, 3, 1, 0.1, 27];
learning_rate = 0.01;

% Define the model function
function dydt_n = model(t, y, params)
    a = params(1);
    b = params(2);
    c = params(3);
    a12_wake = params(4);
    a12_sleep = 2.5 * a12_wake;
    a13_wake = params(5);
    a13_sleep = a * a13_wake;
    a23_sleep = 0.06601;
    a23_wake = a23_sleep/b;
    A_wake = params(6);
    A_sleep = 0.8 * A_wake;
    k_wake = 0.277258;
    k_sleep = c * k_wake;
    
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    dydt_n = zeros(3, 1);
    dydt_n(1) = A_wake * sw_cycle + A_sleep * (1 - sw_cycle) - (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle) + a13_wake * sw_cycle + a13_sleep * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle)) * y(1) - (a23_wake * sw_cycle + a23_sleep * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (a13_wake * sw_cycle + a13_sleep * (1 - sw_cycle)) * y(1) + (a23_wake * sw_cycle + a23_sleep * (1 - sw_cycle)) * y(2) - (k_wake * sw_cycle + k_sleep * (1 - sw_cycle)) * y(3);
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

%Function to perform bootstrapping
function [param_means, param_stds, param_ci] = bootstrap_parameters(csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1, csf_std1, csf_std2, csf_std3, plasma_std1, num_bootstraps, confidence_level)
    n = length(csf_conc_exp1);
    bootstrap_params = zeros(num_bootstraps, 6);

    % Initial guess and bounds (use the same as in your main script)
    initial_guess = [1.1, 1.6, 1.6, 0.3, 0.02, 25];
    lb = [1, 1.5, 1.5, 0.1, 0.01, 20];
    ub = [2, 2, 2, 1, 0.1, 27];

    options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 1000);

    for i = 1:num_bootstraps
        % Resample data with replacement
        indices = randi(n, n, 1);
        bootstrap_csf1 = csf_conc_exp1(indices);
        bootstrap_csf2 = csf_conc_exp2(indices);
        bootstrap_csf3 = csf_conc_exp3(indices);
        bootstrap_plasma1 = plasma_conc_exp1(indices);
        bootstrap_csf_std1 = csf_std1(indices);
        bootstrap_csf_std2 = csf_std2(indices);
        bootstrap_csf_std3 = csf_std3(indices);
        bootstrap_plasma_std1 = plasma_std1(indices);

        % Run optimization for this bootstrap sample
        [bootstrap_params(i,:), ~] = fmincon(@(params) objective_function_fmincon(params, bootstrap_csf1, bootstrap_csf2, bootstrap_csf3, bootstrap_plasma1, bootstrap_csf_std1, bootstrap_csf_std2, bootstrap_csf_std3, bootstrap_plasma_std1), ...
                                             initial_guess, [], [], [], [], lb, ub, [], options);
        
        fprintf('Bootstrap iteration %d/%d completed\n', i, num_bootstraps);
    end

    param_means = mean(bootstrap_params);
    param_stds = std(bootstrap_params);

    %confidence interval
    alpha = 1 - confidence_level;
    lower_percentile = alpha/2 * 100;
    upper_percentile = (1 - alpha/2) * 100;
    param_ci = prctile(bootstrap_params, [lower_percentile, upper_percentile]);

    param_names = {'a', 'b', 'c', 'a12_wake', 'a13_wake', 'A_wake'};
    figure;
    for i = 1:6
        subplot(2,3,i);
        histogram(bootstrap_params(:,i));
        title(param_names{i});
        xlabel('Parameter Values');
        ylabel('Frequency');
        hold on;
        plot([param_ci(1,i), param_ci(1,i)], ylim, 'r--');
        plot([param_ci(2,i), param_ci(2,i)], ylim, 'r--');
        hold off;
    end
end

% Defining Objective Function
function [total_error] = objective_function_fmincon(params, exp_csf1, exp_csf2, exp_csf3, exp_plasma1, csf_std1, csf_std2, csf_std3, plasma_std1)
    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
    
    % Extract last 36 hours of data
    csf_last_36hours_data = sol(233600:237200,2);
    plasma_last_36hours_data = sol(233600:237200,3);
    
    % Select data points corresponding to experimental time points (every 2 hours)
    time_indices = 1:200:3601; % This will give 19 points over 36 hours (0 to 36, step 2)
    csf_model = csf_last_36hours_data(time_indices);
    plasma_model = plasma_last_36hours_data(time_indices);
    
    % Calculate average of first 7 points for model data
    avg_C_model = mean(csf_model(1:7));
    avg_P_model = mean(plasma_model(1:7));
    
    % Normalize the model data
    normalised_C_model = (csf_model/avg_C_model)*100;
    normalised_P_model = (plasma_model/avg_P_model)*100;
    
    csf_weights1 = 1./(csf_std1.^2);
    csf_weights2 = 1./(csf_std2.^2);
    csf_weights3 = 1./(csf_std3.^2);
    plasma_weight1 = 1./(plasma_std1.^2);

    % Calculate error for each dataset
    error_C1 = sqrt(sum(csf_weights1 .* (normalised_C_model - exp_csf1).^2) / sum(csf_weights1));
    error_C2 = sqrt(sum(csf_weights2 .* (normalised_C_model - exp_csf2).^2) / sum(csf_weights2));
    error_C3 = sqrt(sum(csf_weights3 .* (normalised_C_model - exp_csf3).^2) / sum(csf_weights3));
    error_P1 = sqrt(sum(plasma_weight1 .* (normalised_P_model - exp_plasma1).^2) / sum(plasma_weight1));
    
    % Total error is the sum of errors from all datasets
    total_error = error_C1 + error_C2 + error_C3 + error_P1;
end

% Set up optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 2000, 'MaxFunctionEvaluations', 10000);

% Define the optimization problem
problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function_fmincon(params, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1, csf_std1, csf_std2, csf_std3, plasma_std1), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

% Run the optimization
[optimized_params, fval] = fmincon(problem);

num_bootstraps = 100;
confidence_level = 0.95;
[param_means, param_stds, param_ci] = bootstrap_parameters(csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1, csf_std1, csf_std2, csf_std3, plasma_std1, num_bootstraps, confidence_level);

% Run the model for 100 days with optimized parameters
[t_100days, sol_100days] = euler(@(t,y) model(t,y,optimized_params), [0, 24*100], [0,600,15.5], 0.01);

% Save the simulation results using high precision
fileID = fopen('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output.csv', 'w');
fprintf(fileID, 'Time,Brain,CSF,Plasma\n');
for i = 1:length(t_100days)
    fprintf(fileID, '%.2f,%.15f,%.15f,%.15f\n', t_100days(i), sol_100days(i, 1), sol_100days(i, 2), sol_100days(i,3));
end
fclose(fileID);

% Print optimized parameters
fprintf('\nOptimized Parameters for Three Datasets:\n');
fprintf('a: %f\n', optimized_params(1));
fprintf('b: %f\n', optimized_params(2));
fprintf('c: %f\n', optimized_params(3));
fprintf('a12_wake: %f\n', optimized_params(4));
fprintf('a13_wake: %f\n', optimized_params(5));
fprintf('A_wake: %f\n', optimized_params(6));

% Print results
param_names = {'a', 'b', 'c', 'a12_wake', 'a13_wake', 'A_wake'};
for i = 1:length(param_names)
    fprintf('%s: mean = %.4f, std = %.4f, 95%% CI = [%.4f, %.4f]\n', ...
        param_names{i}, param_means(i), param_stds(i), param_ci(1,i), param_ci(2,i));
end

% Run model with optimized parameters
[t, y] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 600, 15.5], 0.01);

% Load saved results for comparison
saved_data = readmatrix('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output.csv');
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

cdata_12hours = sol_100days(236100:237300, 2);
cdata_36hours = sol_100days(236100:239700, 2);
pdata_12hours = sol_100days(236100:237300, 3);
pdata_36hours = sol_100days(236100:239700, 3);
csf_mean = mean(cdata_12hours);
plasma_mean = mean(pdata_12hours);
normalise_C = (cdata_36hours/ csf_mean)*100;
normalise_P = (pdata_36hours / plasma_mean)*100;
disp(length(normalise_C));
disp(length(normalise_P));

cdata_12hours_saved = y_saved(236100:237300, 2);
cdata_36hours_saved = y_saved(236100:238700, 2);
pdata_12hours_saved = y_saved(236100:237300, 3);
pdata_36hours_saved = y_saved(236100:239700, 3);
csf_mean_saved = mean(cdata_12hours);
plasma_mean_saved = mean(pdata_12hours);
normalise_C_saved = (cdata_36hours/ csf_mean_saved)*100;
normalise_P_saved = (pdata_36hours / plasma_mean_saved)*100;
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
plot(time_exp, csf1_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
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
plot(time_exp, csf2_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Lucey2022 - CSF');
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
plot(time_exp, csf3_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Liu2018 - CSF');
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
plot(time_exp, plasma1_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Liu Data - Plasma');
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
