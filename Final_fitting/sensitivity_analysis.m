sigma_A = 0.8;
A = 13;
rbc = 1.5;
rcp = 0.0056;
sigma_bc = 2.5;
rbp = (rbc*(1-133*rcp))/(133*rcp);
rp = 0.28;
sigma_p = 2.89;
sigma_bp_values = linspace(0.1, 10, 100);
sigma_cp_values = linspace(0.1, 10, 100);

[SIGMA_BP, SIGMA_CP] = meshgrid(sigma_bp_values, sigma_cp_values);

% B*
sensitivity_B_sigma_bp = -sigma_A * A * rbp ./ (sigma_bp_values * rbp + sigma_bc * rbc).^2;
sensitivity_B_sigma_cp = zeros(size(sigma_cp_values));

% C*
sensitivity_C_sigma_bp = -sigma_bc * rbc * sigma_A * A * rbp ./ (rcp * SIGMA_CP .* (sigma_bc * rbc + SIGMA_BP * rbp).^2);
sensitivity_C_sigma_cp = -sigma_bc * rbc * sigma_A * A ./ (rcp * SIGMA_CP.^2 .* (sigma_bc * rbc + sigma_bp_values * rbp));

% P*
sensitivity_P_sigma_bp = sigma_A * A * rbp ./ ((sigma_bp_values * rbp + sigma_bc * rbc) * sigma_p * rp).^2;
sensitivity_P_sigma_cp = zeros(size(sigma_cp_values));

figure();
set(gcf, 'Position', [100 100 1000 600]);
subplot(2,3,1)
plot(sigma_bp_values, sensitivity_B_sigma_bp, 'LineWidth', 2.5, 'Color', 'red');
title('\partial \rho_{b}^{s*} / \partial \sigma_{bp}', 'FontSize', 14);
grid on;
xlabel('σ_{bp}');
ylabel('Sensitivity');

subplot(2,3,4)
plot(sigma_cp_values, sensitivity_B_sigma_cp, 'LineWidth', 2.5, 'Color', 'red');
title('\partial \rho_{b}^{s*} / \partial \sigma_{cp}', 'FontSize', 14);
grid on;
xlabel('σ_{cp}');
ylabel('Sensitivity');

subplot(2,3,2)
plot(sigma_bp_values, sensitivity_C_sigma_bp(:,1), 'LineWidth', 2.5, 'Color', 'magenta');
title('\partial \rho_{c}^{s*} / \partial \sigma_{bp}', 'FontSize', 14);
grid on;
xlabel('σ_{bp}');
ylabel('Sensitivity');

subplot(2,3,5)
plot(sigma_cp_values, sensitivity_C_sigma_cp(:,1), 'LineWidth', 2.5, 'Color', 'magenta');
title('\partial \rho_{c}^{s*}/ \partial \sigma_{cp}', 'FontSize', 14);
grid on;
xlabel('σ_{cp}');
ylabel('Sensitivity');

subplot(2, 3, 3)
plot(sigma_bp_values, sensitivity_P_sigma_bp, 'LineWidth', 2.5, 'Color', 'blue');
title('\partial \rho_{p}^{s*} / \partial \sigma_{bp}', 'FontSize', 14);
grid on;
xlabel('\sigma_{bp}');
ylabel('Sensitivity');

subplot(2, 3, 6)
plot(sigma_cp_values, sensitivity_P_sigma_cp, 'LineWidth', 2.5, 'Color', 'blue');
title('\partial \rho_{p}^{s*} / \partial \sigma_{cp}', 'FontSize', 14);
grid on;
xlabel('\sigma_{cp}');
ylabel('Sensitivity');

% colour gradient
colormap(hot);
