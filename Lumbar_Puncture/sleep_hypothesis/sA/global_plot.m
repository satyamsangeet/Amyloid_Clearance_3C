function extract_specific_timepoints(input_file, output_file, start_time, end_time, step)
    % Extract data at specific timepoints and save to a new CSV file
    %
    % Parameters:
    %   input_file  - Path to the input CSV file
    %   output_file - Path to save the output CSV file
    %   start_time  - Starting time point (e.g., 2336)
    %   end_time    - Ending time point (e.g., 2372)
    %   step        - Step size between time points (e.g., 2)
    
    try
        % Read the CSV file
        data = readtable(input_file);
        
        % Create a list of target timepoints
        target_times = start_time:step:end_time;
        
        % Check if 'Time' column exists
        if ~ismember('Time', data.Properties.VariableNames)
            warning('Time column not found in %s', input_file);
            return;
        end
        
        % Extract rows for each target time
        extracted_rows = [];
        for i = 1:length(target_times)
            target_time = target_times(i);
            
            % Find the closest time in the dataframe to the target time
            [~, idx] = min(abs(data.Time - target_time));
            
            % Append the closest row
            if isempty(extracted_rows)
                extracted_rows = data(idx, :);
            else
                extracted_rows = [extracted_rows; data(idx, :)];
            end
        end
        
        % Round to 3 decimal places
        extracted_rows.Time = round(extracted_rows.Time, 3);
        extracted_rows.Compartment_1 = round(extracted_rows.Compartment_1, 3);
        extracted_rows.Compartment_2 = round(extracted_rows.Compartment_2, 3);
        extracted_rows.Compartment_3 = round(extracted_rows.Compartment_3, 3);
        
        % Save to a new CSV file
        writetable(extracted_rows, output_file);
        
        fprintf('Successfully extracted %d timepoints from %s to %s\n', ...
                height(extracted_rows), input_file, output_file);
        
    catch ME
        warning('Error processing %s: %s', input_file, ME.message);
    end
end

% Updated model function to optimize only sigma_bp (a), sigma_cp (b), and rbc (a12_wake)
function dydt_n = model1(t, y)
    r_bc = 0.038;
    r_bp = 0.014;
    r_cp = 0.00537;
    sigma_bc = 1.131;
    sigma_bp = 1.768;
    sigma_cp = 6.100;
    sigma_p = 4.253;
    A = 16.203;
    sigma_A = 0.772;
    r_p = 0.427;

    % Switch
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model2(t, y, params)
    A = 16.203;
    r_p = 0.427;
    r_bp = 0.014;
    r_bc = 0.038;
    r_cp = 0.00537;
    sigma_bc = 1.131;
    sigma_bp = 1.768;
    sigma_cp = 6.100;
    sigma_p = 4.253;
    
    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
    	sigma_A = 0.899;
    else
    	sigma_A = 0.772;
    end
    
    % Switch
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

% Defining Euler-Maruyama Method
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

% Run simulations for both models for 100 days
[t_100days_global1, sol_100days_global1] = euler(@(t,y) model1(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_global2, sol_100days_global2] = euler(@(t,y) model2(t,y), [0, 24*100], [0,600,15.5], 0.01);

% Extract data for all compartments from t=2336 to t=2372
% Create time vector for the extraction period
t_extraction = t_100days_global1(233600:237200);

% Extract data for all 3 compartments for both models during the time period of interest
data_model1_comp1 = sol_100days_global1(233600:237200, 1); % Compartment 1, model 1
data_model1_comp2 = sol_100days_global1(233600:237200, 2); % Compartment 2, model 1
data_model1_comp3 = sol_100days_global1(233600:237200, 3); % Compartment 3, model 1

data_model2_comp1 = sol_100days_global2(233600:237200, 1); % Compartment 1, model 2
data_model2_comp2 = sol_100days_global2(233600:237200, 2); % Compartment 2, model 2
data_model2_comp3 = sol_100days_global2(233600:237200, 3); % Compartment 3, model 2

t_extraction_rounded = round(t_extraction, 3);
data_model1_comp1_rounded = round(data_model1_comp1, 3);
data_model1_comp2_rounded = round(data_model1_comp2, 3);
data_model1_comp3_rounded = round(data_model1_comp3, 3);
data_model2_comp1_rounded = round(data_model2_comp1, 3);
data_model2_comp2_rounded = round(data_model2_comp2, 3);
data_model2_comp3_rounded = round(data_model2_comp3, 3);

% Create tables for each model with the exact column structure requested
model1_table = table(t_extraction_rounded, data_model1_comp1_rounded, data_model1_comp2_rounded, data_model1_comp3_rounded, ...
                    'VariableNames', {'Time', 'Compartment_1', 'Compartment_2', 'Compartment_3'});
model2_table = table(t_extraction_rounded, data_model2_comp1_rounded, data_model2_comp2_rounded, data_model2_comp3_rounded, ...
                    'VariableNames', {'Time', 'Compartment_1', 'Compartment_2', 'Compartment_3'});

% Write tables to CSV files
writetable(model1_table, 'global_model1_h3.csv');
writetable(model2_table, 'global_model2_h3.csv');

% ============================================================================
% PREPROCESSING: Extract specific timepoints for AIC calculation
% ============================================================================

% Define preprocessing parameters
start_time = 2336;
end_time = 2372;
step = 2;

% Extract timepoints for Model 1
fprintf('\nExtracting timepoints for Model 1...\n');
extract_specific_timepoints('global_model1_h3.csv', ...
                            fullfile('global_model1_baseline.csv'), ...
                            start_time, end_time, step);

% Extract timepoints for Model 2
fprintf('Extracting timepoints for Model 2...\n');
extract_specific_timepoints('global_model2_h3.csv', ...
                            fullfile('global_model2_sA.csv'), ...
                            start_time, end_time, step);

% ============================================================================

% For backward compatibility, keep existing cdata_36hours extraction
cdata_36hours1_global1 = sol_100days_global1(233600:237200, 2);
cdata_36hours1_global2 = sol_100days_global2(233600:237200, 2);
pdata_36hours1_global1 = sol_100days_global1(233600:237200, 3);
pdata_36hours1_global2 = sol_100days_global2(233600:237200, 3);

csf_data_file1 = 'data/blattner2020_csf_concentration.csv';
csf_data_file2 = 'data/lucey2018_csf_concentration.csv';
csf_data_file3 = 'data/liu2022_csf_concentration.csv';
plasma_data_file1 = 'data/liu2022_plasma_concentration.csv';
csf_data1 = readtable(csf_data_file1);
csf_data2 = readtable(csf_data_file2);
csf_data3 = readtable(csf_data_file3);
plasma_data1 = readtable(plasma_data_file1);

% Extract data from both plasma files
time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Concentration;
csf_conc_exp2 = csf_data2.Concentration;
csf_conc_exp3 = csf_data3.Concentration;
plasma_conc_exp1 = plasma_data1.Concentration;
csf_lsd1 = csf_data1.LSD;
csf_lsd2 = csf_data2.LSD;
csf_lsd3 = csf_data3.LSD;
plasma_lsd1 = plasma_data1.LSD;
csf_usd1 = csf_data1.USD;
csf_usd2 = csf_data2.USD;
csf_usd3 = csf_data3.USD;
plasma_usd1 = plasma_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;
csf_std2 = (csf_usd2 - csf_lsd2)/2;
csf_std3 = (csf_usd3 - csf_lsd3)/2;
plasma_std1 = (plasma_usd1 - plasma_lsd1)/2;
exp_csf1 = csf_conc_exp1(:);
exp_csf2 = csf_conc_exp2(:);
exp_csf3 = csf_conc_exp3(:);
exp_plasma1 = plasma_conc_exp1(:);

time_indices = 1:200:3601;

% Function to calculate NRMSE
function err = calculate_wrmse(model_data, exp_data)
    errors = (model_data - exp_data).^2;
    error = sqrt(sum(errors) / length(exp_data));
    err = error/abs(max(exp_data) - min(exp_data));
end

csf_36hr_model_global1 = cdata_36hours1_global1(time_indices);
csf_36hr_model_global2 = cdata_36hours1_global2(time_indices);
plasma_36hr_model_global1 = pdata_36hours1_global1(time_indices);
plasma_36hr_model_global2 = pdata_36hours1_global2(time_indices);

wrmse1_global1 = calculate_wrmse(csf_36hr_model_global1, exp_csf1);
wrmse1_global2 = calculate_wrmse(csf_36hr_model_global1, exp_csf2);
wrmse1_global3 = calculate_wrmse(csf_36hr_model_global1, exp_csf3);
wrmse1_global4 = calculate_wrmse(plasma_36hr_model_global1, exp_plasma1);

wrmse1_global5 = calculate_wrmse(csf_36hr_model_global2, exp_csf1);
wrmse1_global6 = calculate_wrmse(csf_36hr_model_global2, exp_csf2);
wrmse1_global7 = calculate_wrmse(csf_36hr_model_global2, exp_csf3);
wrmse1_global8 = calculate_wrmse(plasma_36hr_model_global2, exp_plasma1);

new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);
colormap_jet = jet(5);

x = 1:3601;

figure();
x1 = [time_exp1; flipud(time_exp1)];
inBetween1 = [csf_lsd1; flipud(csf_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Blattner Global (NRMSE: %.3f)', wrmse1_global1));
plot(x, cdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Blattner Global H (NRMSE: %.3f)', wrmse1_global5));
plot(time_exp1, csf_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global1 = interp1(x, cdata_36hours1_global1, time_exp1);
model_points_global2 = interp1(x, cdata_36hours1_global2, time_exp1);
errorbar(time_exp1, csf_conc_exp1, csf_std1, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1, csf_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2372]);
xticks(0:200:3600);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 3600]);
hold off;

figure();
x1 = [time_exp1; flipud(time_exp1)];
inBetween1 = [csf_lsd2; flipud(csf_usd2)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Lucey Global (NRMSE: %.3f)', wrmse1_global2));
plot(x, cdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Lucey Global H (NRMSE: %.3f)', wrmse1_global6));
plot(time_exp1, csf_conc_exp2, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global1 = interp1(x, cdata_36hours1_global1, time_exp1);
model_points_global2 = interp1(x, cdata_36hours1_global2, time_exp1);
errorbar(time_exp1, csf_conc_exp2, csf_std2, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1, csf_conc_exp2, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2372]);
xticks(0:200:3600);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 3600]);
hold off;

figure();
x1 = [time_exp1; flipud(time_exp1)];
inBetween1 = [csf_lsd3; flipud(csf_usd3)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Liu CSF Global (NRMSE: %.3f)', wrmse1_global3));
plot(x, cdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Liu CSF Global H (NRMSE: %.3f)', wrmse1_global7));
plot(time_exp1, csf_conc_exp3, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global1 = interp1(x, cdata_36hours1_global1, time_exp1);
model_points_global2 = interp1(x, cdata_36hours1_global2, time_exp1);
errorbar(time_exp1, csf_conc_exp3, csf_std3, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1, csf_conc_exp3, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2372]);
xticks(0:200:3600);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 3600]);
hold off;

figure();
x1 = [time_exp1; flipud(time_exp1)];
inBetween1 = [plasma_lsd1; flipud(plasma_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, pdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Liu Plasma Global (NRMSE: %.3f)', wrmse1_global4));
plot(x, pdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Liu Plasma Global H (NRMSE: %.3f)', wrmse1_global8));
plot(time_exp1, plasma_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global1 = interp1(x, pdata_36hours1_global1, time_exp1);
model_points_global2 = interp1(x, pdata_36hours1_global2, time_exp1);
errorbar(time_exp1, plasma_conc_exp1, plasma_std1, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1, plasma_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2372]);
xticks(0:200:3600);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 3600]);
hold off;

% Create additional plots to visualize all compartments for both models
figure();
subplot(3,1,1);
plot(t_extraction, data_model1_comp1, 'LineWidth', 2.0, 'Color', 'b', 'DisplayName', 'Model 1');
hold on;
plot(t_extraction, data_model2_comp1, 'LineWidth', 2.0, 'Color', 'r', 'DisplayName', 'Model 2');
title('Compartment 1');
xlabel('Time (hr)');
ylabel('Concentration');
legend('show');
hold off;

subplot(3,1,2);
plot(t_extraction, data_model1_comp2, 'LineWidth', 2.0, 'Color', 'b', 'DisplayName', 'Model 1');
hold on;
plot(t_extraction, data_model2_comp2, 'LineWidth', 2.0, 'Color', 'r', 'DisplayName', 'Model 2');
title('Compartment 2');
xlabel('Time (hr)');
ylabel('Concentration');
legend('show');
hold off;

subplot(3,1,3);
plot(t_extraction, data_model1_comp3, 'LineWidth', 2.0, 'Color', 'b', 'DisplayName', 'Model 1');
hold on;
plot(t_extraction, data_model2_comp3, 'LineWidth', 2.0, 'Color', 'r', 'DisplayName', 'Model 2');
title('Compartment 3');
xlabel('Time (hr)');
ylabel('Concentration');
legend('show');
hold off;

% Print confirmation message
fprintf('Data for model1 and model2 from t=2336 to t=2372 has been extracted and saved\n');
fprintf('CSV files saved: model1_data_t2336_t2372.csv, model2_data_t2336_t2372.csv\n');

f1 = (wrmse1_global1+wrmse1_global2+wrmse1_global3+wrmse1_global4)/4;
f2 = (wrmse1_global5+wrmse1_global6+wrmse1_global7+wrmse1_global8)/4;
final = ((f2 - f1)/f1)*100;
disp(round(final));

function results = calculate_combined_aic_aicc_bic( ...
    y_obs_comp21, y_obs_comp22, y_obs_comp23, y_obs_comp3, ...
    y_pred_comp2, y_pred_comp3, ...
    num_params)

% Calculate combined AIC, AICc, and BIC for two compartments
%
% Inputs:
%   y_obs_comp2  - Observed CSF values
%   y_obs_comp3  - Observed Plasma values
%   y_pred_comp2 - Predicted CSF values
%   y_pred_comp3 - Predicted Plasma values
%   num_params   - Number of model parameters (excluding variance)
%
% Output:
%   results - struct containing AIC/AICc/BIC and diagnostics

    % Ensure column vectors
    y_obs_comp21  = y_obs_comp21(:);
    y_obs_comp22  = y_obs_comp22(:);
    y_obs_comp23  = y_obs_comp23(:);
    y_obs_comp3  = y_obs_comp3(:);
    y_pred_comp2 = y_pred_comp2(:);
    y_pred_comp3 = y_pred_comp3(:);

    % Residuals
    residuals_comp21 = y_obs_comp21 - y_pred_comp2;
    residuals_comp22 = y_obs_comp22 - y_pred_comp2;
    residuals_comp23 = y_obs_comp23 - y_pred_comp2;
    residuals_comp3 = y_obs_comp3 - y_pred_comp3;

    % Combine residuals
    %combined_residuals = [residuals_comp2; residuals_comp3];

    % RSS
    rss1 = sum(residuals_comp21.^2);
    rss2 = sum(residuals_comp22.^2);
    rss3 = sum(residuals_comp23.^2);
    rss4 = sum(residuals_comp3.^2);
    %rss = sum(combined_residuals.^2);

    % Number of observations
    n1 = length(y_obs_comp21);
    n2 = length(y_obs_comp22);
    n3 = length(y_obs_comp23);
    n4 = length(y_obs_comp3);
    n = n1+n2+n3+n4;

    % Number of parameters (+1 for variance)
    k = num_params + 4;

    % MLE variance
    sigma2_mle1 = rss1 / n1;
    sigma2_mle2 = rss2 / n2;
    sigma2_mle3 = rss3 / n3;
    sigma2_mle4 = rss4 / n4;

    % Log-likelihood (Gaussian)
    log_likelihood1 = -n1/2 * (log(2*pi) + log(sigma2_mle1) + 1);
    log_likelihood2 = -n2/2 * (log(2*pi) + log(sigma2_mle2) + 1);
    log_likelihood3 = -n3/2 * (log(2*pi) + log(sigma2_mle3) + 1);
    log_likelihood4 = -n4/2 * (log(2*pi) + log(sigma2_mle4) + 1);
    log_likelihood = log_likelihood1+log_likelihood2+log_likelihood3+log_likelihood4;

    % AIC
    aic_total = 2*k - 2*log_likelihood;

    % AICc
    if (n - k - 1) > 0
        aicc_total = aic_total + (2*k*(k+1)) / (n - k - 1);
    else
        aicc_total = Inf;
    end

    % BIC
    bic_total = k*log(n) - 2*log_likelihood;

    % Per-observation metrics (diagnostic only)
    aic_per_obs  = aic_total  / n;
    aicc_per_obs = aicc_total / n;
    bic_per_obs  = bic_total  / n;

    % Output struct
    results = struct( ...
        'aic_total', aic_total, ...
        'aicc_total', aicc_total, ...
        'bic_total', bic_total, ...
        'aic_per_obs', aic_per_obs, ...
        'aicc_per_obs', aicc_per_obs, ...
        'bic_per_obs', bic_per_obs, ...
        'n_obs', n, ...
        'log_likelihood', log_likelihood ...
    );
end

y_obs_csf1    = csf_conc_exp1;
y_obs_csf2    = csf_conc_exp2;
y_obs_csf3    = csf_conc_exp3;
y_obs_plasma = plasma_conc_exp1;

model_file = 'global_model2_sA.csv';

% Read the preprocessed model data
model_data = readtable(model_file);

y_pred_csf    = model_data.Compartment_2;
y_pred_plasma = model_data.Compartment_3;

num_params = 1;

results = calculate_combined_aic_aicc_bic( ...
    y_obs_csf1, y_obs_csf2, y_obs_csf3, y_obs_plasma, ...
    y_pred_csf, y_pred_plasma, ...
    num_params);

disp(results)