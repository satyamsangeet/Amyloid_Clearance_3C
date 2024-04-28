% Read experimental data from CSV files
data_exp1 = csvread('/home/satyam/Documents/huang_csf42.csv', 1, 0);
data_exp2 = csvread('/home/satyam/Documents/huang_plasma42.csv', 1, 0);
time_exp = data_exp1(:, 1);
exp_csf = data_exp1(:, 2);
exp_plasma = data_exp2(:, 2);

% Define the model function
function dydt = model(t, y, a21, a12_wake, a12_sleep, A_wake, A_sleep, k)
    dydt = zeros(2, 1);
    C = y(1);
    B = y(2);

    % Defining the sleep-wake switch
    sw_cycle = (mod(t, 24) >= 8) & (mod(t, 24) < 24);
    dydt(1) = A_wake * sw_cycle + A_sleep * (1 - sw_cycle) + a21 * y(2) - (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle)) * y(1);
    dydt(2) = (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle)) * y(1) - (a21 + k) * y(2);
end

% Define objective function
function total_error = objective_function(params, exp_csf, exp_plasma)
    a21 = params(1);
    a12_wake = params(2);
    a12_sleep = params(3);
    A_wake = params(4);
    A_sleep = params(5);
    k = params(6);
    
    opts = odeset('MaxStep',0.001);
    [t, sol] = ode45(@(t, y) model(t, y, a21, a12_wake, a12_sleep, A_wake, A_sleep, k), [0, 24*10], [exp_csf(1), exp_plasma(1)], opts);

    % Check if sol has enough data points
    if size(sol, 1) < 36001
        error('Integration did not proceed far enough to generate sufficient data points.');
    end
    
    csf_last_36hours_data = sol(end-36000:end, 1);
    plasma_last_36hours_data = sol(end-36000:end, 2);

    time_indices = 1:36;
    csf_last_36hours = csf_last_36hours_data(time_indices);
    plasma_last_36hours = plasma_last_36hours_data(time_indices);

    avg_C = mean(csf_last_36hours(1:12));
    avg_B = mean(plasma_last_36hours(1:12));

    baseline_C = (csf_last_36hours - avg_C) / avg_C;
    baseline_B = (plasma_last_36hours - avg_B) / avg_B;

    error_C = sum((baseline_C - exp_csf).^2);
    error_B = sum((baseline_B - exp_plasma).^2);
    total_error = error_C + error_B;
end

% Define initial guess for parameters
initial_guess = [0.1, 0.1, 0.3, 30, 1, 2];

% Define bounds for parameters
lb = [0, 0, 0, 0, 0, 0];
ub = [100, 100, 100, 100, 100, 100];

% Minimize the objective function to find optimal parameters using pattern search
options = optimoptions(@patternsearch,'MaxFunctionEvaluations',1000,'MaxIterations',1000);
[result, ~] = patternsearch(@(params) objective_function(params, exp_csf, exp_plasma), initial_guess, [], [], [], [], lb, ub, options);

% Get optimized parameters
a21_opt = result(1);
a12_wake_opt = result(2);
a12_sleep_opt = result(3);
A_wake_opt = result(4);
A_sleep_opt = result(5);
k_opt = result(6);

% Print optimized parameters
fprintf('Optimized Parameters:\n');
fprintf('A_wake: %f\n', A_wake_opt);
fprintf('A_sleep: %f\n', A_sleep_opt);
fprintf('a12_wake: %f\n', a12_wake_opt);
fprintf('a12_sleep: %f\n', a12_sleep_opt);
fprintf('a21: %f\n', a21_opt);
fprintf('k: %f\n', k_opt);

% Solve the model with optimized parameters at each hour and save the results
fileID = fopen('model_output_each_hour_with_mse_loss_bound100_patternsearch.csv', 'w');
fprintf(fileID, 'Time,CSF,Plasma\n');
for hour = 0:120
    try
        opts = odeset('MaxStep',0.001);  % Define opts inside the loop
        [t, sol] = ode45(@(t, y) model(t, y, a21_opt, a12_wake_opt, a12_sleep_opt, A_wake_opt, A_sleep_opt, k_opt), [0 hour], [exp_csf(1), exp_plasma(1)], opts);
        fprintf(fileID, '%d,%f,%f\n', hour, sol(end, 1), sol(end, 2));
    catch e
        fprintf('Error occurred at hour %d: %s\n', hour, e.message);
    end
end
fclose(fileID);