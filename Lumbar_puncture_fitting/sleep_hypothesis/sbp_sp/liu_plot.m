function dydt_n = model1(t, y)
    A = 12.330;
    r_p = 0.337;
    sigma_A = 0.787;
    r_bc = 2.361;
    r_cp = 0.00458;
    r_bp = 0.655;
    sigma_bc = 2.519;
    sigma_cp = 5.540;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        sigma_bp = 6.152;
        sigma_p = 4.154;
    else
        sigma_bp = 6.152;
        sigma_p = 4.154;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model2(t, y)
    A = 12.330;
    r_p = 0.337;
    sigma_A = 0.787;
    r_bc = 2.361;
    r_cp = 0.00458;
    r_bp = 0.655;
    sigma_bc = 2.519;
    sigma_cp = 5.540;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        sigma_bp = 2.467;
        sigma_p = 3.597;
    else
        sigma_bp = 6.152;
        sigma_p = 4.154;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    % Define the ODE system
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

% Save as CSV files with the specified structure: Time, Compartment 1, Compartment 2, Compartment 3
% Create tables for each model with the exact column structure requested
model1_table = table(t_extraction, data_model1_comp1, data_model1_comp2, data_model1_comp3, ...
                    'VariableNames', {'Time', 'Compartment_1', 'Compartment_2', 'Compartment_3'});
model2_table = table(t_extraction, data_model2_comp1, data_model2_comp2, data_model2_comp3, ...
                    'VariableNames', {'Time', 'Compartment_1', 'Compartment_2', 'Compartment_3'});

% Write tables to CSV files
writetable(model1_table, 'liu_model1_h9.csv');
writetable(model2_table, 'liu_model2_h9.csv');

cdata_36hours1_global1 = sol_100days_global1(233600:237200, 2);
cdata_36hours1_global2 = sol_100days_global2(233600:237200, 2);
pdata_36hours1_global1 = sol_100days_global1(233600:237200, 3);
pdata_36hours1_global2 = sol_100days_global2(233600:237200, 3);

csf_data_file1 = 'data/liu2022_csf_concentration.csv';
plasma_data_file1 = 'data/liu2022_plasma_concentration.csv';
csf_data1 = readtable(csf_data_file1);
plasma_data1 = readtable(plasma_data_file1);

% Extract data from both plasma files
time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Concentration;
plasma_conc_exp1 = plasma_data1.Concentration;
csf_lsd1 = csf_data1.LSD;
csf_usd1 = csf_data1.USD;
plasma_lsd1 = plasma_data1.LSD;
plasma_usd1 = plasma_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;
plasma_std1 = (plasma_usd1 - plasma_lsd1)/2;
exp_csf1 = csf_conc_exp1(:);
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
wrmse1_global2 = calculate_wrmse(csf_36hr_model_global2, exp_csf1);
wrmse1_global3 = calculate_wrmse(plasma_36hr_model_global1, exp_plasma1);
wrmse1_global4 = calculate_wrmse(plasma_36hr_model_global2, exp_plasma1);

new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);
colormap_jet = jet(5);

x = 1:3601;

figure();
x1 = [time_exp1/10; flipud(time_exp1/10)];
inBetween1 = [csf_lsd1; flipud(csf_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('H1 (NRMSE: %.3f)', wrmse1_global1));
plot(x, cdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('H1 (NRMSE: %.3f)', wrmse1_global2));
plot(time_exp1/10, csf_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global1 = interp1(x, cdata_36hours1_global1, time_exp1/10);
model_points_global2 = interp1(x, cdata_36hours1_global2, time_exp1/10);
errorbar(time_exp1/10, csf_conc_exp1, csf_std1, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1/10, csf_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1/10, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1/10, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
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
x1 = [time_exp1/10; flipud(time_exp1/10)];
inBetween1 = [plasma_lsd1; flipud(plasma_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, pdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Global Fit (NRMSE: %.3f)', wrmse1_global3));
plot(x, pdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Individual Fit (NRMSE: %.3f)', wrmse1_global4));
model_points_global3 = interp1(x, pdata_36hours1_global1, time_exp1/10);
model_points_global4 = interp1(x, pdata_36hours1_global2, time_exp1/10);
errorbar(time_exp1/10, plasma_conc_exp1, plasma_std1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot(time_exp1/10, plasma_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
scatter(time_exp1/10, plasma_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1/10, model_points_global3, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp1/10, model_points_global4, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Individual Fit Points');
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