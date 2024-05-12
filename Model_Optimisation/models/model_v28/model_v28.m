% Read experimental data from CSV files
data_exp1 = csvread('/home/satyam/Documents/huang_csf42.csv', 1, 0);
data_exp2 = csvread('/home/satyam/Documents/huang_plasma42.csv', 1, 0);
time_exp = data_exp1(:, 1);
exp_csf = data_exp1(:, 2);
exp_plasma = data_exp2(:, 2);

% Define global variable for loss_values
global loss_values;
loss_values = [];

function dydt_n = model(t, y, a12_wake, a13_wake, a23_wake, A_wake, A_sleep, a, k)
    sw_cycle = (mod(t,24) >= 8 && mod(t,24) < 24);
    dydt_n = zeros(3, 1);
    dydt_n(1) = A_wake * sw_cycle + (A_sleep * (1 - sw_cycle)) - (a13_wake * sw_cycle + a*a13_wake*(1-sw_cycle) + a12_wake*sw_cycle + a*a12_wake*(1 - sw_cycle))*y(1);
    dydt_n(2) = (a12_wake*sw_cycle + a*a12_wake*(1-sw_cycle))*y(1) - (a23_wake*sw_cycle + a*a23_wake*(1-sw_cycle))*y(2);
    dydt_n(3) = (a23_wake*sw_cycle + a*a23_wake*(1-sw_cycle))*y(2) + (a13_wake*sw_cycle + a*a13_wake*(1-sw_cycle))*y(1) - k*y(3);
end

% Define objective function
function total_error = objective_function(params, exp_csf, exp_plasma)
    global loss_values;

    a12_wake = params(1);
    a13_wake = params(2);
    a23_wake = params(3);
    A_wake = params(4);
    A_sleep = params(5);
    a = params(6);
    k = params(7);
    
    opts = odeset('MaxStep',0.01);
    [t, sol] = ode45(@(t, y) model(t, y, a12_wake, a13_wake, a23_wake, A_wake, A_sleep, a, k), [0, 24*10], [exp_csf(1), exp_plasma(1), 0], opts);
    
    csf_last_36hours_data = sol(20000:23600, 1);
    plasma_last_36hours_data = sol(20000:23600, 2);

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
initial_guess = [0.1, 0.1, 0.2, 10, 7, 0.9, 0.2];

% Define bounds for parameters
lb = [0, 0, 0, 10, 7, 0.89, 0.1];
ub = [0.1, 0.1, 0.3, 10.5, 7.5, 0.99, 0.6];

% Minimize the objective function to find optimal parameters using pattern search
options = optimoptions(@patternsearch,'MaxFunctionEvaluations',300,'MaxIterations',300);
[result, ~] = patternsearch(@(params) objective_function(params, exp_csf, exp_plasma), initial_guess, [], [], [], [], lb, ub, options);

% Get optimized parameters
a12_wake_opt = result(1);
a13_wake_opt = result(2);
a23_wake_opt = result(3);
A_wake_opt = result(4);
A_sleep_opt = result(5);
a_opt = result(6);
k_opt = result(7);

% Print optimized parameters
fprintf('Optimized Parameters:\n');
fprintf('A_wake: %f\n', A_wake_opt);
fprintf('A_sleep: %f\n', A_sleep_opt);
fprintf('a12_wake: %f\n', a12_wake_opt);
fprintf('k: %f\n', k_opt);
fprintf('a: %f\n', a_opt);
fprintf('a13_wake: %f\n', a13_wake_opt);
fprintf('a23_wake: %f\n', a23_wake_opt);

%Saving the loss values
lossID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v28_loss.csv', 'w');
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

% Solve the model with optimized parameters at each hour and save the results
fileID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v28_output.csv', 'w');
fprintf(fileID, 'Time,CSF,Plasma\n');
for hour = 0:0.1:240
    try
        opts = odeset('MaxStep',0.01);  % Define opts inside the loop
        tspan = [0, hour];  % Define tspan for each iteration
        [t, sol] = ode45(@(t, y) model(t, y, a12_wake_opt, a13_wake_opt, a23_wake_opt, A_wake_opt, A_sleep_opt, a_opt, k_opt), tspan, [exp_csf(1), exp_plasma(1),0], opts);
        fprintf(fileID, '%d,%f,%f\n', hour, sol(end, 1), sol(end, 2));
    catch e
        fprintf('Error occurred at hour %d: %s\n', hour, e.message);
    end
end
fclose(fileID);