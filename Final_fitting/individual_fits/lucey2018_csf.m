csf_data_file2 = '/home/satyam/Documents/combined_fitting_data/final_fitting/lucey2018_csf42.csv';
csf_data2 = readtable(csf_data_file2);
time_exp2 = csf_data2.Time;
csf_conc_exp2 = csf_data2.Conc;
csf_lsd2 = csf_data2.LSD;
csf_usd2 = csf_data2.USD;
csf_std2 = (csf_usd2 - csf_lsd2)/2;

initial_guess = [0.005, 3, 7, 7];  
lb = [0, 1, 1, 1]; 
ub = [0.0074, 10, 10, 10];
learning_rate = 0.01;

function dydt_n = model(t, y, params)
    r_cp = params(1);
    sigma_bp = params(2);
    sigma_cp = params(3);
    sigma_p = params(4);

    A = 13;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_p = 0.28;
    r_bp = (r_bc*(1-133*r_cp))/133*r_cp;

    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_cp * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

% Euler-Maruyama Method
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

function [param_means, param_stds, param_ci] = bootstrap_parameters(csf_conc_exp2, csf_std2, num_bootstraps, confidence_level)
    n = length(csf_conc_exp2);
    bootstrap_params = zeros(num_bootstraps, 4);

    initial_guess = [0.005, 3, 7, 7];  
    lb = [0, 1, 1, 1]; 
    ub = [0.0074, 10, 10, 10];

    options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 1000);

    for i = 1:num_bootstraps
        indices = randi(n, n, 1);
        bootstrap_csf2 = csf_conc_exp2(indices);
        bootstrap_csf_std2 = csf_std2(indices);

        [bootstrap_params(i,:), ~] = fmincon(@(params) objective_function_fmincon(params, bootstrap_csf2, bootstrap_csf_std2), ...
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

    param_names = {'r_cp', 'sigma_bp', 'sigma_cp', 'sigma_p'};
    figure;
    for i = 1:4
        subplot(2,2,i);
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

function [total_error] = objective_function_fmincon(params, exp_csf2, csf_std2)
    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
    
    csf_last_36hours_data = sol(233600:237200,2);
    time_indices = 1:200:3601; % This will give 19 points over 36 hours (0 to 36, step 2)
    csf_model = csf_last_36hours_data(time_indices);
    avg_C_model = mean(csf_model(1:7));
    normalised_C_model = (csf_model/avg_C_model)*100;
    csf_weights2 = 1./(csf_std2.^2);
    error_C2 = sqrt(sum(csf_weights2 .* (normalised_C_model - exp_csf2).^2) / sum(csf_weights2));
    total_error = error_C2;
end

global optimization_history;
optimization_history = [];

function stop = save_loss_values(x, optimValues, state)
    global optimization_history;
    stop = false;
    
    if strcmp(state, 'iter')
        optimization_history = [optimization_history; optimValues.iteration, optimValues.fval];
    end
end

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 10000, ...
    'OutputFcn', @save_loss_values, ...
    'Algorithm','interior-point', ...
    'FiniteDifferenceStepSize',1e-6, ...
    'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-6);

problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function_fmincon(params, csf_conc_exp1, csf_std1), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

num_starts = 10;
md = MultiStart('UseParallel', true, 'Display', 'iter');

% Run the optimization
[optimized_params, fval] = run(md, problem, num_starts);

num_bootstraps = 100;
confidence_level = 0.95;
[param_means, param_stds, param_ci] = bootstrap_parameters(csf_conc_exp2, csf_std2, num_bootstraps, confidence_level);

[t_100days, sol_100days] = euler(@(t,y) model(t,y,optimized_params), [0, 24*100], [0,600,15.5], 0.01);

fileID = fopen('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values_lucey2018.csv', 'w');
fprintf(fileID, 'Time,Brain,CSF,Plasma\n');
for i = 1:length(t_100days)
    fprintf(fileID, '%.2f,%.15f,%.15f,%.15f\n', t_100days(i), sol_100days(i, 1), sol_100days(i, 2), sol_100days(i,3));
end
fclose(fileID);

global optimization_history;
writematrix(optimization_history, 'optimization_loss_history.csv', 'Delimiter', ',');

fprintf('\nOptimized Parameters for Three Datasets:\n');
fprintf('r_cp: %f\n', optimized_params(1));
fprintf('sigma_bp: %f\n', optimized_params(2));
fprintf('sigma_cp: %f\n', optimized_params(3));
fprintf('sigma_p: %f\n', optimized_params(3));

param_names = {'r_bc', 'sigma_bp', 'sigma_cp', 'sigma_p'};
for i = 1:length(param_names)
    fprintf('%s: mean = %.4f, std = %.4f, 95%% CI = [%.4f, %.4f]\n', ...
        param_names{i}, param_means(i), param_stds(i), param_ci(1,i), param_ci(2,i));
end

[t, y] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 600, 15.5], 0.01);

saved_data = readmatrix('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values_lucey2018.csv');
t_saved = saved_data(:, 1);
y_saved = saved_data(:, 2:end);

figure;
plot(optimization_history(:,1), optimization_history(:,2), 'b-', 'LineWidth', 2.0);
hold on;
scatter(optimization_history(:,1), optimization_history(:,2), 50, 'blue', 'filled');
xlabel('Iteration');
ylabel('Loss values');
grid on;

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
csf_mean = mean(cdata_12hours);
normalise_C = (cdata_36hours/ csf_mean)*100;
disp(length(normalise_C));

cdata_12hours_saved = y_saved(233600:234800, 2);
cdata_36hours_saved = y_saved(233600:237200, 2);
csf_mean_saved = mean(cdata_12hours_saved);
normalise_C_saved = (cdata_36hours_saved / csf_mean_saved)*100;
disp(length(normalise_C_saved));

experimental_csf2 = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/lucey2018_csf42.csv');
experimental_csf2.Time = experimental_csf2.Time - 100;
time_exp = experimental_csf2.Time;
csf2_conc_exp = experimental_csf2.Conc;

x = 1:3601;
ticks = 101:100:3601;
x = 1:3601;
disp(length(x));

new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);
valid_indices = ticks + 1 <= length(normalise_C);
valid_ticks = ticks(valid_indices);

figure;
subplot(3, 1, 1);
hold on;
plot(x, normalise_C, 'r', 'LineWidth', 2.0, 'HandleVisibility','off');
plot(x, normalise_C_saved, 'b--', 'DisplayName', 'Fitted Data');
plot(time_exp, csf2_conc_exp, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Lucey2018 - CSF');
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
