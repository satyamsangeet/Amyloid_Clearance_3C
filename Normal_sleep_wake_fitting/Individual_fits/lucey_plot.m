% global fit
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

%individual fit
function dydt_n = model2(t, y)
    r_bc = 0.062;
    r_bp = 0.040;
    r_cp = 0.0156;
    sigma_bc = 1.002;
    sigma_bp = 2.479;
    sigma_cp = 4.087;
    sigma_p = 4;
    A = 55.780;
    sigma_A = 0.485;
    r_p = 0.300;

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

[t_100days_global1, sol_100days_global1] = euler(@(t,y) model1(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_global2, sol_100days_global2] = euler(@(t,y) model2(t,y), [0, 24*100], [0,600,15.5], 0.01);

cdata_36hours1_global1 = sol_100days_global1(233600:237600, 2);
cdata_36hours1_global2 = sol_100days_global2(233600:237600, 2);

csf_data_file1 = 'data/lucey_wake_conc.csv';
csf_data1 = readtable(csf_data_file1);

% Extract data from both plasma files
time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Concentration;
csf_lsd1 = csf_data1.LSD;
csf_usd1 = csf_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;

exp_csf1 = csf_conc_exp1(:);

time_indices = 1:200:4001;
selected_indices1 = [1:9, 13:21];

% Function to calculate NRMSE
function err = calculate_wrmse(model_data, exp_data)
    disp(size(exp_data));
    disp(size(model_data));
    errors = (model_data - exp_data).^2;
    error = sqrt(sum(errors) / length(exp_data));
    err = error/abs(max(exp_data) - min(exp_data));
end

csf_36hr_model_global1 = cdata_36hours1_global1(time_indices);
norm_csf_global11 = csf_36hr_model_global1(selected_indices1);
wrmse1_global1 = calculate_wrmse(norm_csf_global11, exp_csf1);

csf_36hr_model_global2 = cdata_36hours1_global2(time_indices);
norm_csf_global12 = csf_36hr_model_global2(selected_indices1);
wrmse1_global2 = calculate_wrmse(norm_csf_global12, exp_csf1);

time_start = 2330;
time_end = 2380;
time_interval = 2;

selected_time_indices = find(mod(t_100days_global1, time_interval) == 0 & ...
                            t_100days_global1 >= time_start & ...
                            t_100days_global1 <= time_end);

time_selected = t_100days_global1(selected_time_indices);
c1_selected = sol_100days_global1(selected_time_indices, 1);
c2_selected = sol_100days_global1(selected_time_indices, 2);
c3_selected = sol_100days_global1(selected_time_indices, 3);

simulation_data = [time_selected, c1_selected, c2_selected, c3_selected];
csv_filename = 'model_simulation_data.csv';
writematrix(simulation_data, csv_filename);
disp(['Simulation data saved to ', csv_filename]);

ticks = 101:100:4001;
x = 1:4001;

disp(size(x));
disp(size(cdata_36hours1_global1));

% Define new ticks and labels
new_ticks = 1:100:4000+100;
new_labels = 1:length(new_ticks);
colormap_jet = jet(2);

% Plot
figure();
x1 = [time_exp1; flipud(time_exp1)];
inBetween1 = [csf_lsd1; flipud(csf_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('Global Fit (NRMSE: %.3f)', wrmse1_global1));
plot(x, cdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('Individual Fit (NRMSE: %.3f)', wrmse1_global2));
model_points_global1 = interp1(x, cdata_36hours1_global1, time_exp1);
model_points_global2 = interp1(x, cdata_36hours1_global2, time_exp1);
errorbar(time_exp1, csf_conc_exp1, csf_std1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'lucey2020 - CSF', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
plot(time_exp1, csf_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'lucey2020 - CSF');
scatter(time_exp1, csf_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Global Fit Points');
scatter(time_exp1, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Individual Fit Points');
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


% Residual analysis
fprintf('time_exp1 range: %.2f to %.2f\n', min(time_exp1), max(time_exp1))
fprintf('x range: %.2f to %.2f\n', min(x), max(x))
fprintf('time_exp1 first few values: ')
disp(time_exp1(1:5)')
fprintf('x first few values: ')
disp(x(1:5))

time_exp1_adjusted = time_exp1 + 1;
model_csf_lucey = interp1(x, cdata_36hours1_global2, time_exp1_adjusted);
fprintf('NaNs in lucey: %d\n', sum(isnan(model_csf_lucey)))

residuals_lucey = exp_csf1 - model_csf_lucey;
time_hours_plot = time_exp1 / 100;

fprintf('Lucey CSF  - Mean: %.3f, SD: %.3f pg/ml\n', mean(residuals_lucey), std(residuals_lucey))

% Plot
figure();
stem(time_hours_plot, residuals_lucey, 'filled', ...
    'Color', [0.2 0.4 0.8], 'LineWidth', 1.5)
hold on
yline(0, 'k--', 'LineWidth', 1.5)
xlabel('Time (hr)', 'FontWeight', 'bold')
ylabel('Residual (pg/ml)', 'FontWeight', 'bold')
title('lucey et al. 2020 - CSF', 'FontWeight', 'bold')
text(1, max(abs(residuals_lucey))*0.85, ...
    sprintf('Mean = %.2f pg/ml', mean(residuals_lucey)), 'FontSize', 9)
grid on
hold off
