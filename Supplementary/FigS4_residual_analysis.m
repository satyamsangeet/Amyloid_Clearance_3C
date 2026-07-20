function dydt_n = model_blattner(t, y)
    r_bc = 0.019;
    r_bp = 0.034;
    r_cp = 0.01540;
    sigma_bc = 1.660;
    sigma_bp = 1.816;
    sigma_cp = 5.740;
    sigma_p = 3.610;
    A = 84.523;
    sigma_A = 0.633;
    r_p = 0.298;
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model_lucey(t, y)
    r_bc = 0.062;
    r_bp = 0.040;
    r_cp = 0.01560;
    sigma_bc = 1.002;
    sigma_bp = 2.479;
    sigma_cp = 4.087;
    sigma_p = 4.000;
    A = 55.780;
    sigma_A = 0.485;
    r_p = 0.300;
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model_liu(t, y)
    r_bc = 0.015;
    r_bp = 0.014;
    r_cp = 0.00320;
    sigma_bc = 1.110;
    sigma_bp = 1.297;
    sigma_cp = 6.552;
    sigma_p = 3.055;
    A = 14.450;
    sigma_A = 0.750;
    r_p = 0.475;
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

[t_ind_b, sol_ind_b] = euler(@(t,y) model_blattner(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_ind_l, sol_ind_l] = euler(@(t,y) model_lucey(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_ind_liu, sol_ind_liu] = euler(@(t,y) model_liu(t,y), [0, 24*100], [0,600,15.5], 0.01);

time_exp1_adjusted = time_exp1 + 1;

model_csf_blattner_ind = interp1(x, sol_ind_b(233600:237600, 2), time_exp1_adjusted);
model_csf_lucey_ind = interp1(x, sol_ind_l(233600:237600, 2), time_exp1_adjusted);
model_csf_liu_ind = interp1(x, sol_ind_liu(233600:237600, 2), time_exp1_adjusted);
model_plasma_liu_ind = interp1(x, sol_ind_liu(233600:237600, 3), time_exp1_adjusted);

res_blattner_ind = exp_csf1 - model_csf_blattner_ind;
res_lucey_ind = exp_csf2 - model_csf_lucey_ind;
res_liu_csf_ind = exp_csf3 - model_csf_liu_ind;
res_liu_plasma_ind = exp_plasma1 - model_plasma_liu_ind;

time_hours_plot = time_exp1 / 100;

fprintf('Blattner CSF  - Mean: %.3f, SD: %.3f pg/ml\n', mean(res_blattner_ind,'omitnan'), std(res_blattner_ind,'omitnan'))
fprintf('Lucey CSF     - Mean: %.3f, SD: %.3f pg/ml\n', mean(res_lucey_ind,'omitnan'), std(res_lucey_ind,'omitnan'))
fprintf('Liu CSF       - Mean: %.3f, SD: %.3f pg/ml\n', mean(res_liu_csf_ind,'omitnan'), std(res_liu_csf_ind,'omitnan'))
fprintf('Liu Plasma    - Mean: %.3f, SD: %.3f pg/ml\n', mean(res_liu_plasma_ind,'omitnan'), std(res_liu_plasma_ind,'omitnan'))

% Plot
figure();
subplot(2,2,1)
stem(time_hours_plot, res_blattner_ind, 'filled', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5)
hold on
yline(0, 'k--', 'LineWidth', 1.5)
xlabel('Time (hr)', 'FontWeight', 'bold')
ylabel('Residual (pg/ml)', 'FontWeight', 'bold')
title('Blattner et al. 2020 - CSF (Individual Fit)', 'FontWeight', 'bold')
text(1, max(abs(res_blattner_ind))*0.85, sprintf('Mean = %.2f pg/ml', mean(res_blattner_ind,'omitnan')), 'FontSize', 9)
grid on
hold off

subplot(2,2,2)
stem(time_hours_plot, res_lucey_ind, 'filled', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5)
hold on
yline(0, 'k--', 'LineWidth', 1.5)
xlabel('Time (hr)', 'FontWeight', 'bold')
ylabel('Residual (pg/ml)', 'FontWeight', 'bold')
title('Lucey et al. 2018 - CSF (Individual Fit)', 'FontWeight', 'bold')
text(1, max(abs(res_lucey_ind))*0.85, sprintf('Mean = %.2f pg/ml', mean(res_lucey_ind,'omitnan')), 'FontSize', 9)
grid on
hold off

subplot(2,2,3)
stem(time_hours_plot, res_liu_csf_ind, 'filled', 'Color', [0.2 0.7 0.3], 'LineWidth', 1.5)
hold on
yline(0, 'k--', 'LineWidth', 1.5)
xlabel('Time (hr)', 'FontWeight', 'bold')
ylabel('Residual (pg/ml)', 'FontWeight', 'bold')
title('Liu et al. 2023 - CSF (Individual Fit)', 'FontWeight', 'bold')
text(1, max(abs(res_liu_csf_ind))*0.85, sprintf('Mean = %.2f pg/ml', mean(res_liu_csf_ind,'omitnan')), 'FontSize', 9)
grid on
hold off

subplot(2,2,4)
stem(time_hours_plot, res_liu_plasma_ind, 'filled', 'Color', [0.8 0.5 0.1], 'LineWidth', 1.5)
hold on
yline(0, 'k--', 'LineWidth', 1.5)
xlabel('Time (hr)', 'FontWeight', 'bold')
ylabel('Residual (pg/ml)', 'FontWeight', 'bold')
title('Liu et al. 2023 - Plasma (Individual Fit)', 'FontWeight', 'bold')
text(1, max(abs(res_liu_plasma_ind))*0.85, sprintf('Mean = %.2f pg/ml', mean(res_liu_plasma_ind,'omitnan')), 'FontSize', 9)
grid on
hold off

sgtitle('Residual Analysis: Individual Fits', 'FontWeight', 'bold', 'FontSize', 12)
