function dydt_n = model1(t, y)
    A = 11.724;
    sigma_bc = 1.185;
    sigma_cp = 4.173;
    sigma_A = 0.731;
    sigma_bp = 6.607;
    sigma_p = 4.028;
    r_bp = 0.242;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        r_bc = 1.319;
        r_cp = 0.00545;
        r_p = 0.317;
    else
        r_bc = 1.319;
        r_cp = 0.00545;
        r_p = 0.317;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model2(t, y)
    A = 11.724;
    r_p = 0.317;
    sigma_bc = 1.185;
    sigma_cp = 4.173;
    sigma_A = 0.731;
    sigma_bp = 6.607;
    sigma_p = 4.028;
    r_bp = 0.242;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        r_bc = 1.630;
        r_cp = 0.00176;
        r_p = 0.119;
    else
        r_bc = 1.319;
        r_cp = 0.00545;
        r_p = 0.317;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model3(t, y)
    A = 11.724;
    r_p = 0.317;
    sigma_bc = 1.185;
    sigma_cp = 4.173;
    sigma_A = 0.731;
    sigma_bp = 6.607;
    sigma_p = 4.028;
    r_bp = 0.242;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        r_bc = 1.329;
        r_cp = 0.00229;
        r_p = 0.161;
    else
        r_bc = 1.319;
        r_cp = 0.00545;
        r_p = 0.317;
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
[t_100days_global3, sol_100days_global3] = euler(@(t,y) model3(t,y), [0, 24*100], [0,600,15.5], 0.01);

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

data_model3_comp1 = sol_100days_global3(233600:237200, 1); % Compartment 1, model 2
data_model3_comp2 = sol_100days_global3(233600:237200, 2); % Compartment 2, model 2
data_model3_comp3 = sol_100days_global3(233600:237200, 3); % Compartment 3, model 2

% Save as CSV files with the specified structure: Time, Compartment 1, Compartment 2, Compartment 3
% Create tables for each model with the exact column structure requested
model1_table = table(t_extraction, data_model1_comp1, data_model1_comp2, data_model1_comp3, ...
                    'VariableNames', {'Time', 'Compartment_1', 'Compartment_2', 'Compartment_3'});
model2_table = table(t_extraction, data_model2_comp1, data_model2_comp2, data_model2_comp3, ...
                    'VariableNames', {'Time', 'Compartment_1', 'Compartment_2', 'Compartment_3'});
model3_table = table(t_extraction, data_model3_comp1, data_model3_comp2, data_model3_comp3, ...
                    'VariableNames', {'Time', 'Compartment_1', 'Compartment_2', 'Compartment_3'});

% Write tables to CSV files
writetable(model1_table, 'global_model.csv');
writetable(model2_table, 'global_error1.csv');
writetable(model3_table, 'global_error2.csv');

cdata_36hours1_global1 = sol_100days_global1(233600:237200, 2);
pdata_36hours1_global1 = sol_100days_global1(233600:237200, 3);
cdata_36hours2_global2 = sol_100days_global2(233600:237200, 2);
pdata_36hours2_global2 = sol_100days_global2(233600:237200, 3);
cdata_36hours3_global3 = sol_100days_global3(233600:237200, 2);
pdata_36hours3_global3 = sol_100days_global3(233600:237200, 3);

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
time_exp2 = csf_data2.Time;
time_exp3 = csf_data3.Time;
csf_conc_exp1 = csf_data1.Concentration;
csf_conc_exp2 = csf_data2.Concentration;
csf_conc_exp3 = csf_data3.Concentration;
plasma_conc_exp1 = plasma_data1.Concentration;
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
plasma_std = (plasma_usd1 - plasma_lsd1)/2;

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
csf_36hr_model_global2 = cdata_36hours2_global2(time_indices);
csf_36hr_model_global3 = cdata_36hours3_global3(time_indices);
plasma_36hr_model_global1 = pdata_36hours1_global1(time_indices);
plasma_36hr_model_global2 = pdata_36hours2_global2(time_indices);
plasma_36hr_model_global3 = pdata_36hours3_global3(time_indices);

wrmse1_global1 = calculate_wrmse(csf_36hr_model_global1, exp_csf1);
wrmse1_global2 = calculate_wrmse(csf_36hr_model_global2, exp_csf1);
wrmse1_global3 = calculate_wrmse(csf_36hr_model_global3, exp_csf1);
wrmse1_global4 = calculate_wrmse(csf_36hr_model_global1, exp_csf2);
wrmse1_global5 = calculate_wrmse(csf_36hr_model_global2, exp_csf2);
wrmse1_global6 = calculate_wrmse(csf_36hr_model_global3, exp_csf2);
wrmse1_global7 = calculate_wrmse(csf_36hr_model_global1, exp_csf3);
wrmse1_global8 = calculate_wrmse(csf_36hr_model_global2, exp_csf3);
wrmse1_global9 = calculate_wrmse(csf_36hr_model_global3, exp_csf3);
wrmse1_global10 = calculate_wrmse(plasma_36hr_model_global1, exp_plasma1);
wrmse1_global11 = calculate_wrmse(plasma_36hr_model_global2, exp_plasma1);
wrmse1_global12 = calculate_wrmse(plasma_36hr_model_global3, exp_plasma1);

ticks = 101:100:3601;
x = 1:3601;

disp(size(x));
disp(size(cdata_36hours1_global1));

% Define new ticks and labels
new_ticks = 1:10036200+100;
new_labels = 1:length(new_ticks);
colormap_jet = jet(3);

% Plot comparison CSF
figure();
x1 = [time_exp1; flipud(time_exp1)];
inBetween1 = [csf_lsd1; flipud(csf_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Global (NRMSE: %.3f)', wrmse1_global1));
plot(x, cdata_36hours2_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Error 1 (NRMSE: %.3f)', wrmse1_global2));
plot(x, cdata_36hours3_global3, 'LineWidth', 2.0, 'Color', colormap_jet(3,:), 'DisplayName', sprintf('Error 2 (NRMSE: %.3f)', wrmse1_global3));
plot(time_exp1, csf_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global11 = interp1(x, cdata_36hours1_global1, time_exp1);
model_points_global12 = interp1(x, cdata_36hours2_global2, time_exp1);
model_points_global13 = interp1(x, cdata_36hours3_global3, time_exp1);
errorbar(time_exp1, csf_conc_exp1, csf_std1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1, csf_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1, model_points_global11, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp1, model_points_global12, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp1, model_points_global13, 80, colormap_jet(3,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2386]);
xticks(0:200:4000);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 4000]);
hold off;

figure();
x1 = [time_exp2; flipud(time_exp2)];
inBetween1 = [csf_lsd2; flipud(csf_usd2)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Global (NRMSE: %.3f)', wrmse1_global4));
plot(x, cdata_36hours2_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Error 1 (NRMSE: %.3f)', wrmse1_global5));
plot(x, cdata_36hours3_global3, 'LineWidth', 2.0, 'Color', colormap_jet(3,:), 'DisplayName', sprintf('Error 2 (NRMSE: %.3f)', wrmse1_global6));
plot(time_exp2, csf_conc_exp2, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global21 = interp1(x, cdata_36hours1_global1, time_exp2);
model_points_global22 = interp1(x, cdata_36hours2_global2, time_exp2);
model_points_global23 = interp1(x, cdata_36hours3_global3, time_exp2);
errorbar(time_exp2, csf_conc_exp2, csf_std2, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp2, csf_conc_exp2, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp2, model_points_global21, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp2, model_points_global22, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp2, model_points_global23, 80, colormap_jet(3,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2386]);
xticks(0:200:4000);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 4000]);
hold off;

figure();
x1 = [time_exp3/10; flipud(time_exp3/10)];
inBetween1 = [csf_lsd3; flipud(csf_usd3)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Global (NRMSE: %.3f)', wrmse1_global7));
plot(x, cdata_36hours2_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Error 1 (NRMSE: %.3f)', wrmse1_global8));
plot(x, cdata_36hours3_global3, 'LineWidth', 2.0, 'Color', colormap_jet(3,:), 'DisplayName', sprintf('Error 2 (NRMSE: %.3f)', wrmse1_global9));
plot(time_exp3/10, csf_conc_exp3, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global31 = interp1(x, cdata_36hours1_global1, time_exp3/10);
model_points_global32 = interp1(x, cdata_36hours2_global2, time_exp3/10);
model_points_global33 = interp1(x, cdata_36hours3_global3, time_exp3/10);
errorbar(time_exp3/10, csf_conc_exp3, csf_std3, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp3/10, csf_conc_exp3, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp3/10, model_points_global31, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp3/10, model_points_global32, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp3/10, model_points_global33, 80, colormap_jet(3,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2386]);
xticks(0:200:4000);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 4000]);
hold off;

figure();
x1 = [time_exp3/10; flipud(time_exp3/10)];
inBetween1 = [plasma_lsd1; flipud(plasma_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, pdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Global (NRMSE: %.3f)', wrmse1_global10));
plot(x, pdata_36hours2_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Error 1 (NRMSE: %.3f)', wrmse1_global11));
plot(x, pdata_36hours3_global3, 'LineWidth', 2.0, 'Color', colormap_jet(3,:), 'DisplayName', sprintf('Error 2 (NRMSE: %.3f)', wrmse1_global12));
plot(time_exp3/10, plasma_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global41 = interp1(x, pdata_36hours1_global1, time_exp3/10);
model_points_global42 = interp1(x, pdata_36hours2_global2, time_exp3/10);
model_points_global43 = interp1(x, pdata_36hours3_global3, time_exp3/10);
errorbar(time_exp3/10, plasma_conc_exp1, plasma_std, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp3/10, plasma_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp3/10, model_points_global41, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp3/10, model_points_global42, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp3/10, model_points_global43, 80, colormap_jet(3,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2386]);
xticks(0:200:4000);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 4000]);
hold off;


% Plot comparison
figure;
subplot(3, 1, 1);
hold on;
plot(t_100days_global1, sol_100days_global1(:, 1), 'r', 'LineWidth', 3.0, 'DisplayName', 'Brain Compartment - Raw');
legend('show');
ylabel('Amyloid Concentration', 'FontSize', 12, 'FontWeight','bold');
xlim([2336, 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2360, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2376, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2384, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
hold off;

subplot(3, 1, 2);
hold on;
plot(t_100days_global1, sol_100days_global1(:, 2), 'r', 'LineWidth', 3.0, 'DisplayName', 'CSF Compartment');
legend('show');
ylabel('Amyloid Concentration', 'FontSize', 12, 'FontWeight','bold');
xlim([2336, 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2360, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2376, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2384, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
hold off;

subplot(3, 1, 3);
hold on;
plot(t_100days_global1, sol_100days_global1(:, 3), 'r', 'LineWidth', 3.0, 'DisplayName', 'Plasma Compartment');
legend('show');
xlabel('Time (hr)', 'FontSize', 12, 'FontWeight','bold');
ylabel('Amyloid Concentration', 'FontSize', 12, 'FontWeight','bold');
xlim([2336, 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2360, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2376, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2384, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
hold off;