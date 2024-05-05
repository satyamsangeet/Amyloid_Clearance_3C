% Read experimental data from CSV files
data_exp1 = csvread('/home/satyam/Documents/huang_csf42.csv', 1, 0);
data_exp2 = csvread('/home/satyam/Documents/huang_plasma42.csv', 1, 0);
time_exp = data_exp1(:, 1);
exp_csf = data_exp1(:, 2);
exp_plasma = data_exp2(:, 2);

% Define the model function
function dydt = model(t, y, a12_wake, a12_sleep, A_wake, A_sleep, k)
    dydt = zeros(2, 1);
    C = y(1);
    B = y(2);

    % Defining the sleep-wake switch
    sw_cycle = (mod(t, 24) >= 8) & (mod(t, 24) < 24);
    dydt(1) = A_wake * sw_cycle + A_sleep * (1 - sw_cycle) - (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle)) * y(1);
    dydt(2) = (a12_wake * sw_cycle + a12_sleep * (1 - sw_cycle)) * y(1) - k * y(2);
end

loss_values = [];

% Define objective function
function total_error = objective_function(params, exp_csf, exp_plasma, loss_values)
    a12_wake = params(1);
    a12_sleep = params(2);
    A_wake = params(3);
    A_sleep = params(4);
    k = params(5);
    
    opts = odeset('MaxStep',0.001);
    [t, sol] = ode45(@(t, y) model(t, y, a12_wake, a12_sleep, A_wake, A_sleep, k), [0, 24*10], [exp_csf(1), exp_plasma(1)], opts);

    % Check if sol has enough data points
    if size(sol, 1) < 204001
        error('Integration did not proceed far enough to generate sufficient data points.');
    end
    
    csf_last_36hours_data = sol(end-204000:end, 1);
    plasma_last_36hours_data = sol(end-204000:end, 2);

    time_indices = 1:36;
    csf_last_36hours = csf_last_36hours_data(time_indices);
    plasma_last_36hours = plasma_last_36hours_data(time_indices);

    avg_C = mean(csf_last_36hours(1:12));
    avg_B = mean(plasma_last_36hours(1:12));

    baseline_C = (csf_last_36hours - avg_C) / avg_C;
    baseline_B = (plasma_last_36hours - avg_B) / avg_B;

    error_C = sqrt(mean((baseline_C - exp_csf).^2));
    error_B = sqrt(mean((baseline_B - exp_plasma).^2));
    total_error = error_C + error_B;
    loss_values(end+1) = total_error;
end

% Define initial guess for parameters
initial_guess = [0.05, 0.3, 15, 0.1, 0.2];

% Define bounds for parameters
lb = [0, 0, 10, 0, 0];
ub = [0.1, 0.5, 20, 5, 1];

% Minimize the objective function to find optimal parameters using pattern search
options = optimoptions(@patternsearch,'MaxFunctionEvaluations',5000,'MaxIterations',5000);
[result, ~] = patternsearch(@(params) objective_function(params, exp_csf, exp_plasma, loss_values), initial_guess, [], [], [], [], lb, ub, options);

% Get optimized parameters
a12_wake_opt = result(1);
a12_sleep_opt = result(2);
A_wake_opt = result(3);
A_sleep_opt = result(4);
k_opt = result(5);

% Print optimized parameters
fprintf('Optimized Parameters:\n');
fprintf('A_wake: %f\n', A_wake_opt);
fprintf('A_sleep: %f\n', A_sleep_opt);
fprintf('a12_wake: %f\n', a12_wake_opt);
fprintf('a12_sleep: %f\n', a12_sleep_opt);
fprintf('k: %f\n', k_opt);

%Saving the loss values
%lossID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v2_loss.csv', 'w');
%fprintf(lossID, 'Iteration,Loss\n');
%for i = 1:length(loss_values)
%    try
%        fprintf(lossID, '%d,%f\n', i, loss_values(i));
%    catch
%        fprintf('Error occurred while saving loss value at iteration %d\n', i);
%    end
%end
%fclose(lossID);

% Solve the model with optimized parameters at each hour and save the results
fileID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v5_output.csv', 'w');
fprintf(fileID, 'Time,CSF,Plasma\n');
for hour = 0:240
    try
        opts = odeset('MaxStep',0.001);  % Define opts inside the loop
        [t, sol] = ode45(@(t, y) model(t, y, a12_wake_opt, a12_sleep_opt, A_wake_opt, A_sleep_opt, k_opt), [0 hour], [exp_csf(1), exp_plasma(1)], opts);
        fprintf(fileID, '%d,%f,%f\n', hour, sol(end, 1), sol(end, 2));
    catch e
        fprintf('Error occurred at hour %d: %s\n', hour, e.message);
    end
end
fclose(fileID);