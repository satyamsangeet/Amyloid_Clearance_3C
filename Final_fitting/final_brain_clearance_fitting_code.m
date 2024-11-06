% Model Fitting

% Data file path
csf_data_file1 = '/home/satyam/Documents/combined_fitting_data/final_fitting/blattner2020_csf42.csv';
csf_data_file2 = '/home/satyam/Documents/combined_fitting_data/final_fitting/lucey2018_csf42.csv';
csf_data_file3 = '/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_csf42.csv';
plasma_data_file1 = '/home/satyam/Documents/combined_fitting_data/final_fitting/liu2022_plasma42.csv';

csf_data1 = readtable(csf_data_file1);
csf_data2 = readtable(csf_data_file2);
csf_data3 = readtable(csf_data_file3);
plasma_data1 = readtable(plasma_data_file1);

% Extracting time, conc, lower bound and upper bound of SD
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

% Initail guess, lower bound, upper bound
initial_guess = [3, 12, 0.003];  % Optimizing sigma_bp (a), sigma_cp (b), rbc (a12_wake)
lb = [1, 10, 0.0001];  % Lower bounds for sigma_bp, sigma_cp, and rbc
ub = [10, 15, 0.005];  % Upper bounds for sigma_bp, sigma_cp, and rbc
learning_rate = 0.01;

% Setting the optim options
options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 2000, 'MaxFunctionEvaluations', 10000);

% Running optim
[optimized_params, fval] = fmincon(@(params) objective_function_fmincon(params, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1, csf_std1, csf_std2, csf_std3, plasma_std1), ...
    initial_guess, [], [], [], [], lb, ub, @nlcon, options);

% Running bootstrapping
num_bootstraps = 10;
confidence_level = 0.95;
[param_means, param_stds, param_ci] = bootstrap_parameters(csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1, csf_std1, csf_std2, csf_std3, plasma_std1, num_bootstraps, confidence_level);

% 2400 hours or 100 days run (24 hours * 100 days)
[t_100days, sol_100days] = euler(@(t,y) model(t,y,optimized_params), [0, 24*100], [0,600,15.5], 0.01);

% Saving data
fileID = fopen('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values.csv', 'w');
fprintf(fileID, 'Time,Brain,CSF,Plasma\n');
for i = 1:length(t_100days)
    fprintf(fileID, '%.2f,%.15f,%.15f,%.15f\n', t_100days(i), sol_100days(i, 1), sol_100days(i, 2), sol_100days(i,3));
end
fclose(fileID);

% Optim params
fprintf('\nOptimized Parameters for Three Datasets:\n');
fprintf('sigma_bp: %f\n', optimized_params(1));
fprintf('sigma_p: %f\n', optimized_params(2));
fprintf('r_cp: %f\n', optimized_params(3));

param_names = {'sigma_bp', 'sigma_p', 'r_cp'};
for i = 1:length(param_names)
    fprintf('%s: mean = %.4f, std = %.4f, 95%% CI = [%.4f, %.4f]\n', ...
        param_names{i}, param_means(i), param_stds(i), param_ci(1,i), param_ci(2,i));
end

% Running model with optim params
[t, y] = euler(@(t, y) model(t, y, optimized_params), [0, 24*100], [0, 600, 15.5], 0.01);

saved_data = readmatrix('/home/satyam/Documents/combined_fitting_data/final_fitting/combined_model_v1_output_new_values.csv');
t_saved = saved_data(:, 1);
y_saved = saved_data(:, 2:end);

createVisualizationPlots(t, y, t_saved, y_saved, sol_100days, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1);

function [c, ceq] = nlcon(params)
    r_cp = params(3);
    sigma_cp = 0.077/r_cp; % Calculate sigma_cp based on r_cp

    % Define bounds for sigma_cp
    lower_bound = 0.0578/r_cp;
    upper_bound = 0.077/r_cp;

    % Constraints: c <= 0
    c = [sigma_cp - upper_bound;
         lower_bound - sigma_cp];
    
    ceq = []; % No equality constraints
end

function dydt_n = model(t, y, params)
    sigma_bp = params(1);
    sigma_p = params(2);
    r_cp = params(3);
    sigma_cp = 0.077/r_cp;

    % Derived parameters
    A = 14;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_bp = (r_bc*(1-155*r_cp))/(155*r_cp);
    r_p = 0.277258;

    % Switch between sleep and wake states
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_cp * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

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

function [param_means, param_stds, param_ci] = bootstrap_parameters(csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1, csf_std1, csf_std2, csf_std3, plasma_std1, num_bootstraps, confidence_level)
    n = length(csf_conc_exp1);
    bootstrap_params = zeros(num_bootstraps, 3);

    % Initial guess and bounds
    initial_guess = [3, 12, 0.003];
    lb = [1, 10, 0.0001];
    ub = [10, 15, 0.005];

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
                                             initial_guess, [], [], [], [], lb, ub, @nlcon, options);
        
        fprintf('Bootstrap iteration %d/%d completed\n', i, num_bootstraps);
    end

    param_means = mean(bootstrap_params);
    param_stds = std(bootstrap_params);

    % Calculate confidence interval
    alpha = 1 - confidence_level;
    lower_percentile = alpha/2 * 100;
    upper_percentile = (1 - alpha/2) * 100;
    param_ci = prctile(bootstrap_params, [lower_percentile, upper_percentile]);

    % Create bootstrap distribution plots
    param_names = {'sigma_bp', 'sigma_p', 'r_cp'};
    figure;
    for i = 1:3
        subplot(1,3,i);
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

function createVisualizationPlots(t, y, t_saved, y_saved, sol_100days, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp1) 
    figure; 
    subplot(3,1,1); 
    plot(t, y(:,2), 'b-', 'LineWidth', 1.5); 
    hold on; plot(t_saved, y_saved(:,2), 'r--', 'LineWidth', 1.5); 
    plot(1:2:36, csf_conc_exp1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); 
    xlabel('Time (hours)'); ylabel('CSF Concentration'); 
    legend('Model - Current', 'Saved Model', 'Experimental CSF Data 1'); 
    title('CSF Concentration Over Time');
    
    subplot(3,1,2);
    plot(t, y(:,2), 'b-', 'LineWidth', 1.5); hold on;
    plot(t_saved, y_saved(:,2), 'r--', 'LineWidth', 1.5);
    plot(1:2:36, csf_conc_exp2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    xlabel('Time (hours)'); ylabel('CSF Concentration');
    legend('Model - Current', 'Saved Model', 'Experimental CSF Data 2');
    title('CSF Concentration Over Time');

    subplot(3,1,3);
    plot(t, y(:,2), 'b-', 'LineWidth', 1.5); hold on;
    plot(t_saved, y_saved(:,2), 'r--', 'LineWidth', 1.5);
    plot(1:2:36, csf_conc_exp3, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    xlabel('Time (hours)'); ylabel('CSF Concentration');
    legend('Model - Current', 'Saved Model', 'Experimental CSF Data 3');
    title('CSF Concentration Over Time');
    
    figure;
    plot(t, y(:,3), 'g-', 'LineWidth', 1.5); hold on;
    plot(t_saved, y_saved(:,3), 'm--', 'LineWidth', 1.5);
    plot(1:2:36, plasma_conc_exp1, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    xlabel('Time (hours)'); ylabel('Plasma Concentration');
    legend('Model - Current', 'Saved Model', 'Experimental Plasma Data');
    title('Plasma Concentration Over Time');
end

