rng(42);
global optimization_history;
optimization_history = [];

csf_data_file1 = 'data/blattner2020_csf_concentration.csv';
csf_data_file2 = 'data/lucey2018_csf_concentration.csv';
csf_data_file3 = 'data/liu2022_csf_concentration.csv';
plasma_data_file = 'data/liu2022_plasma_concentration.csv';
csf_data1 = readtable(csf_data_file1);
csf_data2 = readtable(csf_data_file2);
csf_data3 = readtable(csf_data_file3);
plasma_data = readtable(plasma_data_file);

time_exp1 = csf_data1.Time;
time_exp2 = csf_data2.Time;
time_exp3 = csf_data3.Time;
csf_conc_exp1 = csf_data1.Concentration;
csf_conc_exp2 = csf_data2.Concentration;
csf_conc_exp3 = csf_data3.Concentration;
plasma_conc_exp = plasma_data.Concentration;
csf_lsd1 = csf_data1.LSD;
csf_lsd2 = csf_data2.LSD;
csf_lsd3 = csf_data3.LSD;
csf_usd1 = csf_data1.USD;
csf_usd2 = csf_data2.USD;
csf_usd3 = csf_data3.USD;
plasma_lsd = plasma_data.LSD;
plasma_usd = plasma_data.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;
csf_std2 = (csf_usd2 - csf_lsd2)/2;
csf_std3 = (csf_usd3 - csf_lsd3)/2;
plasma_std = (plasma_usd - plasma_lsd)/2;

function dydt_n = model(t, y, params)
    A = 11.724240389400634;
    sigma_A = 0.731113610034996;
    r_bp = 0.242172221021308;
    sigma_bp = 6.607200387445705;
    r_p = 0.317320504848485;
    sigma_p = 4.028397512330345;

    if(t>=2336 && t<2372)
        r_bc = params(1);
        r_cp = params(2);
        sigma_bc = params(3);
        sigma_cp = params(4);
    else
        r_bc = 1.319956380713120;
        sigma_bc = 1.185932416265309;
        r_cp = 0.005459771155221;
        sigma_cp = 4.173247278028787;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function [t, w] = euler(F, endpoints, initial_conditions, ts)
    if length(endpoints) == 2
        h = ts; %delta_t (seconds)
        total_time = endpoints(2) - endpoints(1);
        num_steps = floor(total_time / h);
        t = linspace(endpoints(1), endpoints(2), num_steps + 1); %Creating time vector
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints;
    end
    w = zeros(num_steps+1, length(initial_conditions));
    w(1,:) = initial_conditions;
    for k = 1:num_steps
        w(k+1,:) = w(k,:) + F(t(k), w(k,:))' * h;
    end
    t = t(:);
end

function [total_error] = objective_function_fmincon(params, exp_csf1, exp_csf2, exp_csf3, exp_plasma)
    global optimization_history;
    optimization_history = [optimization_history; params];

    [t, sol] = euler(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
    
    csf_last_36hours_data = sol(233600:237200,2);
    plasma_last_36hours_data = sol(233600:237200,3);
    time_indices = 1:200:3601;
    csf_model = csf_last_36hours_data(time_indices);
    plasma_model = plasma_last_36hours_data(time_indices);

    csf_model = csf_model(:);
    plasma_model = plasma_model(:);
    exp_csf1 = exp_csf1(:);
    exp_csf2 = exp_csf2(:);
    exp_csf3 = exp_csf3(:);
    exp_plasma = exp_plasma(:);

    errors_C1 = (csf_model - exp_csf1) .^ 2;
    errors_C2 = (csf_model - exp_csf2) .^ 2;
    errors_C3 = (csf_model - exp_csf3) .^ 2;
    errors_P1 = (plasma_model - exp_plasma) .^ 2;

    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));
    error_C2 = sqrt(sum(errors_C2) / length(exp_csf2));
    error_C3 = sqrt(sum(errors_C3) / length(exp_csf3));
    error_P1 = sqrt(sum(errors_P1) / length(exp_plasma));

    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));
    nrmse_C2 = error_C2/abs(max(exp_csf2) - min(exp_csf2));
    nrmse_C3 = error_C3/abs(max(exp_csf3) - min(exp_csf3));
    nrmse_P1 = error_P1/abs(max(exp_plasma) - min(exp_plasma));

    total_error = (nrmse_C1+nrmse_C2+nrmse_C3+nrmse_P1)/4;

    fprintf('Parameters: r_bc=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_cp=%.4f\n', params(1), params(2), params(3), params(4));
    fprintf('Individual Study Errors: CSF1=%.4f, CSF2=%.4f, CSF3=%.4f, Plasma=%.4f\n', error_C1, error_C2, error_C3, error_P1);
    fprintf('Total WRMSE: %.4f\n', total_error);
end

function stop = output_function(x, optimValues, state)
    global optimization_history;
    stop = false;
    
    if strcmp(state, 'iter')
        optimization_history = [optimization_history; x];
    end
end

function plot_parameter_space(solutions, num_starts)
    params_matrix = zeros(num_starts, 4);
    fvals = zeros(num_starts, 1);
    for i = 1:num_starts
        params_matrix(i,:) = solutions(i).params;
        fvals(i) = solutions(i).Fval;
    end
    
    param_names = {'r_bc', 'r_cp', 'sigma_bc', 'sigma_cp'};

    colors = jet(num_starts);
    figure('Position', [100, 100, 1200, 800]);
    plot_handles = cell(num_starts, 1);
    
    for i = 1:4
        subplot(3, 4, i);
        hold on;
        
        for j = 1:num_starts
            h = scatter(params_matrix(j,i), fvals(j), 50, colors(j,:), 'filled');
            
            if i == 1
                plot_handles{j} = h;
            end
        end
        
        xlabel(param_names{i});
        ylabel('Objective Value');
        title(sprintf('Effect of %s on Objective', param_names{i}));
        grid on;
        hold off;
    end

    subplot(3, 4, 11);  % Use position 11 for legend
    hold on;

    for j = 1:num_starts
        scatter(NaN, NaN, 50, colors(j,:), 'filled', 'DisplayName', sprintf('Run %d', j));
    end
    
    legend('Location', 'best');
    axis off;
    sgtitle('Parameter Space Exploration');
    saveas(gcf, 'global_fit/parameter_space.png');
    print('global_fit/parameter_space_highres', '-dpng', '-r300');

    for param_idx = 1:4
        param_name = param_names{param_idx};
        csv_filename = sprintf('global_fit/parameter_space_%s.csv', param_name);
        
        csv_data = zeros(num_starts, 3);
        for run_idx = 1:num_starts
            csv_data(run_idx, 1) = run_idx;  % Run number
            csv_data(run_idx, 2) = params_matrix(run_idx, param_idx);  % Parameter value
            csv_data(run_idx, 3) = fvals(run_idx);  % Objective function value
        end
        
        % Save to CSV
        writematrix(csv_data, csv_filename);
        fprintf('Saved parameter space data for %s to %s\n', param_name, csv_filename);
    end
end


num_starts = 20;
initial_guess = [1.5, 0.005, 2.5, 3];
lb = [0.1, 0.001, 0.1, 0.1];
ub = [3, 3, 4, 7];


figure('Position', [100, 100, 1200, 800]);

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'Algorithm', 'interior-point', ...
    'FiniteDifferenceStepSize', 1e-4, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-6, ...
    'OutputFcn', @output_function);

problem = createOptimProblem('fmincon', ...
    'objective', @(params) objective_function_fmincon(params, csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp), ...
    'x0', initial_guess, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

customStartPoints = zeros(num_starts, length(initial_guess));
customStartPoints(1,:) = initial_guess;

for i = 2:num_starts
    customStartPoints(i,:) = lb + rand(1, length(lb)) .* (ub - lb);
end

solutions = struct('params', cell(1, num_starts), 'Fval', cell(1, num_starts), 'Exitflag', cell(1, num_starts), 'loss_history', cell(1, num_starts));
all_trajectories = cell(num_starts, 1);
global optimization_history;

if ~exist('global_fit', 'dir')
    mkdir('global_fit');
end

best_fval = Inf;
best_params = [];

for i = 1:num_starts
    optimization_history = [];
    
    problem.x0 = customStartPoints(i,:);
    [x, fval, exitflag, output] = fmincon(problem);
    all_trajectories{i} = optimization_history;
    
    solutions(i).params = x(:)';  % Store optimized parameters in 'params'
    solutions(i).Fval = fval;
    solutions(i).Exitflag = exitflag;
    
    if fval < best_fval
        best_fval = fval;
        best_params = x;
    end
    
    run_results = struct('params', x(:)', ...
        'fval', fval, ...
        'exitflag', exitflag, ...
        'trajectory', optimization_history);
    save(sprintf('global_fit/run_%d_results.mat', i), 'run_results');
    
    fprintf('\nRun %d:\n', i);
    fprintf('Parameters: r_bc=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_cp=%.4f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3), solutions(i).params(4));
    fprintf('Objective value: %.6f\n', fval);
    fprintf('Completed run %d/%d\n', i, num_starts);
end

all_solutions = zeros(num_starts, length(initial_guess));
all_fvals = zeros(num_starts, 1);
for i = 1:num_starts
    all_solutions(i,:) = solutions(i).params;
    all_fvals(i) = solutions(i).Fval;
end

param_names = {'r_bc', 'r_cp', 'sigma_bc', 'sigma_cp'};
correlation_matrix = corrcoef(all_solutions);

fprintf('\nParameter Correlation Matrix:\n');
fprintf('%-10s', 'Param');
for i = 1:length(param_names)
    fprintf('%-10s', param_names{i});
end
fprintf('\n');

for i = 1:length(param_names)
    fprintf('%-10s', param_names{i});
    for j = 1:length(param_names)
        fprintf('%-10.4f', correlation_matrix(i,j));
    end
    fprintf('\n');
end

fileID_corr = fopen('global_fit/correlation_matrix.txt', 'w');
fprintf(fileID_corr, 'Parameter Correlation Matrix:\n');
fprintf(fileID_corr, '%-10s', 'Param');
for i = 1:length(param_names)
    fprintf(fileID_corr, '%-10s', param_names{i});
end
fprintf(fileID_corr, '\n');

for i = 1:length(param_names)
    fprintf(fileID_corr, '%-10s', param_names{i});
    for j = 1:length(param_names)
        fprintf(fileID_corr, '%-10.4f', correlation_matrix(i,j));
    end
    fprintf(fileID_corr, '\n');
end
fclose(fileID_corr);

%heatmap
figure('Position', [100, 100, 800, 600]);
imagesc(correlation_matrix);
colorbar;
colormap(jet);
title('Parameter Correlation Matrix');
xticks(1:length(param_names));
yticks(1:length(param_names));
xticklabels(param_names);
yticklabels(param_names);
xtickangle(45);
textStrings = num2str(correlation_matrix(:), '%.2f');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:length(param_names));
hText = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));
textColors = repmat(correlation_matrix(:) > midValue, 1, 3);
set(hText, {'Color'}, num2cell(textColors, 2));
saveas(gcf, 'global_fit/correlation_heatmap.png');

final_results = struct('best_params', best_params, ...
    'best_fval', best_fval, ...
    'all_solutions', all_solutions, ...
    'all_fvals', all_fvals, ...
    'all_trajectories', all_trajectories, ...
    'correlation_matrix', correlation_matrix);
save('global_fit/final_results.mat', 'final_results');

fileID = fopen('global_fit/detailed_results.txt', 'w');
fprintf(fileID, 'Optimization Results Summary\n');
fprintf(fileID, '===========================\n\n');

for i = 1:num_starts
    fprintf(fileID, 'Run %d:\n', i);
    fprintf(fileID, 'Parameters: r_bc=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_cp=%.4f\n', ...
        solutions(i).params(1), solutions(i).params(2), solutions(i).params(3), solutions(i).params(4));
    
    [t_run, sol_run] = euler(@(t,y) model(t,y,solutions(i).params), [0, 24*100], [0,600,15.5], 0.01);
    
    csf_last_36hours = sol_run(233600:237200,2);
    plasma_last_36hours = sol_run(233600:237200,3);
    
    time_indices = 1:200:3601;
    csf_model = csf_last_36hours(time_indices);
    plasma_model = plasma_last_36hours(time_indices);
    
    csf_model1 = csf_model(:);
    plasma_model1 = plasma_model(:);
    exp_csf1 = csf_conc_exp1(:);
    exp_csf2 = csf_conc_exp2(:);
    exp_csf3 = csf_conc_exp3(:);
    exp_plasma = plasma_conc_exp(:);

    errors_C1 = (csf_model1 - exp_csf1) .^ 2;
    errors_C2 = (csf_model1 - exp_csf2) .^ 2;
    errors_C3 = (csf_model1 - exp_csf3) .^ 2;
    errors_P1 = (plasma_model1 - exp_plasma) .^ 2;

    error_C1 = sqrt(sum(errors_C1) / length(exp_csf1));
    error_C2 = sqrt(sum(errors_C2) / length(exp_csf2));
    error_C3 = sqrt(sum(errors_C3) / length(exp_csf3));
    error_P1 = sqrt(sum(errors_P1) / length(exp_plasma));

    nrmse_C1 = error_C1/abs(max(exp_csf1) - min(exp_csf1));
    nrmse_C2 = error_C2/abs(max(exp_csf2) - min(exp_csf2));
    nrmse_C3 = error_C3/abs(max(exp_csf3) - min(exp_csf3));
    nrmse_P1 = error_P1/abs(max(exp_plasma) - min(exp_plasma));

    total_error = (nrmse_C1+nrmse_C2+nrmse_C3+nrmse_P1)/4;

    fprintf(fileID, 'Individual Errors:\n');
    fprintf(fileID, '  CSF1 Error: %.6f\n', error_C1);
    fprintf(fileID, '  CSF2 Error: %.6f\n', error_C2);
    fprintf(fileID, '  CSF3 Error: %.6f\n', error_C3);
    fprintf(fileID, '  Plasma Error: %.6f\n', error_P1);
    fprintf(fileID, 'Total Error: %.6f\n', total_error);
end

fprintf(fileID, '\nBest Solution:\n');
fprintf(fileID, 'Parameters: r_bc=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_cp=%.4f\n', best_params(1), best_params(2), best_params(3), best_params(4));
fprintf(fileID, 'Objective value: %.6f\n', best_fval);
fclose(fileID);

for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = arrayfun(@(idx) objective_function_fmincon(loss_curve(idx,:), csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp), iterations);
    loss_data = [iterations(:), loss_values(:)];
    writematrix(loss_data, sprintf('global_fit/loss_curve_run_%d.csv', i));
end

figure;
hold on;
for i = 1:num_starts
    loss_curve = all_trajectories{i};
    iterations = 1:size(loss_curve, 1);
    loss_values = arrayfun(@(idx) objective_function_fmincon(loss_curve(idx,:), csf_conc_exp1, csf_conc_exp2, csf_conc_exp3, plasma_conc_exp), iterations);
    plot(iterations, loss_values, 'DisplayName', sprintf('Run %d', i));
end
xlabel('Iteration');
ylabel('Loss Value');
title('Loss Curves for Each Run');
legend;
hold off;
saveas(gcf, 'global_fit/loss_curves.png');

plot_parameter_space(solutions, num_starts);

fprintf('\nBest solution found:\n');
fprintf('Parameters: r_bc=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_cp=%.4f\n', ...
    best_params(1), best_params(2), best_params(3), best_params(4));
fprintf('Objective value: %.6f\n', best_fval);
