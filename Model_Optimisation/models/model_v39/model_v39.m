% Read experimental data from CSV files
data_exp1 = csvread('/home/satyam/Downloads/Two_compartment/studies/lucey2018/lucey2018_sd.csv', 2, 0);
time_exp = data_exp1(:, 1);
exp_csf = data_exp1(:, 3);

% Define global variable for loss_values
global loss_values;
loss_values = [];

function dydt = model(t, y, a12_wake, a13_wake, a23_wake, A_wake, A_sleep, a, k)
    % Determine sleep-wake cycle based on time t
    if t < 2352
        sw_cycle = (mod(t, 24) >= 8) & (mod(t, 24) < 24);
    elseif t >= 2352 && t < 2376
        sw_cycle = 1;
    else
        adjusted_time = t - 2376;
        sw_cycle = (mod(adjusted_time, 24) >= 8) && (mod(adjusted_time, 24) < 24);
    end

    dydt = zeros(3, 1);
    dydt(1) = A_wake * sw_cycle + (A_sleep * (1 - sw_cycle)) - (a13_wake * sw_cycle + a * a13_wake * (1 - sw_cycle) + a12_wake * sw_cycle + a * a12_wake * (1 - sw_cycle)) * y(1);
    dydt(2) = (a12_wake * sw_cycle + a * a12_wake * (1 - sw_cycle)) * y(1) - (a23_wake * sw_cycle + a * a23_wake * (1 - sw_cycle)) * y(2);
    dydt(3) = (a23_wake * sw_cycle + a * a23_wake * (1 - sw_cycle)) * y(2) + (a13_wake * sw_cycle + a * a13_wake * (1 - sw_cycle)) * y(1) - k * y(3);
end


% Define objective function
function total_error = objective_function(params, exp_csf)
    global loss_values;

    a12_wake = params(1);
    a13_wake = params(2);
    a23_wake = params(3);
    A_wake = params(4);
    A_sleep = params(5);
    a = params(6);
    k = params(7);
    
    opts = odeset('MaxStep',0.01);
    [t, sol] = ode45(@(t, y) model(t, y, a12_wake, a13_wake, a23_wake, A_wake, A_sleep, a, k), [0, 24*10], [0, 0, 0], opts);
    
    csf_last_36hours_data = sol(23360:23720, 2);

    time_indices = 1:2:36;
    csf_last_36hours = csf_last_36hours_data(time_indices);

    avg_C = mean(csf_last_36hours(1:12));

    normalised_C = (csf_last_36hours - avg_C) / avg_C;

    error_C = sqrt(mean((normalised_C - exp_csf).^2));
    total_error = error_C;
    
    % Update the global loss_values for this iteration
    loss_values(end+1) = total_error;
end

% Define initial guess for parameters
initial_guess = [0.325, 0.25, 0.25, 12.5, 8.5, 1, 0.6];

% Define bounds for parameters
lb = [0.1, 0.1, 0.1, 10, 7, 1, 0.3];
ub = [0.5, 0.5, 0.5, 13, 9, 2, 0.8];

% Minimize the objective function to find optimal parameters using pattern search
options = optimoptions(@patternsearch,'MaxFunctionEvaluations',1000,'MaxIterations',1000);
[result, ~] = patternsearch(@(params) objective_function(params, exp_csf), initial_guess, [], [], [], [], lb, ub, options);

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
lossID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v39_loss.csv', 'w');
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
fileID = fopen('/home/satyam/Documents/MATLAB/Two_Compartment_Model/Optimisation/model_v39_output.csv', 'w');
fprintf(fileID, 'Time,Brain,CSF,Plasma\n');
for hour = 0:0.1:240
    try
        opts = odeset('MaxStep',0.01);  % Define opts inside the loop
        tspan = [0, hour];  % Define tspan for each iteration
        [t, sol] = ode45(@(t, y) model(t, y, a12_wake_opt, a13_wake_opt, a23_wake_opt, A_wake_opt, A_sleep_opt, a_opt, k_opt), tspan, [0, 0, 0], opts);
        fprintf(fileID, '%d,%f,%f,%f\n', hour, sol(end, 1), sol(end, 2), sol(end, 3));
    catch e
        fprintf('Error occurred at hour %d: %s\n', hour, e.message);
    end
end
fclose(fileID);