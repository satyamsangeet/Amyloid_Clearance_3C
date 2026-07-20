set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 0.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 8, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 8, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);

set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

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

% Defining Euler-Maruyama Method
function [t, w] = euler(F, endpoints, initial_conditions, ts)
    % Calculate number of steps and time vector
    total_time = endpoints(2) - endpoints(1);
    num_steps = floor(total_time / ts);
    
    % Initialize time vector
    t = linspace(endpoints(1), endpoints(2), num_steps + 1)';
    
    % Initialize solution matrix
    w = zeros(length(t), length(initial_conditions));
    w(1,:) = initial_conditions;
    
    % Euler iteration
    for k = 1:num_steps
        w(k+1,:) = w(k,:) + F(t(k), w(k,:))' * ts;
    end
end

% simulation
[t_100days_gwrmse1, sol_100days_gwrmse1] = euler(@(t,y) model1(t,y), [0, 24*100], [0,600,20], 0.01);

% Parameters for clearance calculations
r_bc = 0.038;
r_bp = 0.014;
r_cp = 0.00537;
sigma_bc = 1.131;
sigma_bp = 1.768;
sigma_cp = 6.100;

analysis_idx = find(t_100days_gwrmse1 >= 2336 & t_100days_gwrmse1 <= 2360);
t_analysis = t_100days_gwrmse1(analysis_idx);
sol_analysis = sol_100days_gwrmse1(analysis_idx, :);
time_of_day = mod(t_analysis, 24);
wake_idx = (time_of_day >= 8) & (time_of_day < 24);
sleep_idx = (time_of_day >= 0) & (time_of_day < 8);

% mean concentrations
wake_brain_conc = mean(sol_analysis(wake_idx, 1));
sleep_brain_conc = mean(sol_analysis(sleep_idx, 1));
wake_csf_conc = mean(sol_analysis(wake_idx, 2));
sleep_csf_conc = mean(sol_analysis(sleep_idx, 2));

fprintf('STEADY-STATE BRAIN CONCENTRATIONS:\n');
fprintf('  During wake:  %.2f pg/mL\n', wake_brain_conc);
fprintf('  During sleep: %.2f pg/mL\n\n', sleep_brain_conc);

fprintf('WAKE PERIOD:\n');
fprintf('  Brain to CSF efflux (r_bc):           %.4f h⁻¹\n', r_bc);
fprintf('  Blood-Brain Barrier (r_bp):           %.4f h⁻¹\n', r_bp);
fprintf('  Blood-CSF Barrier (r_cp):           %.4f h⁻¹\n', r_cp);
fprintf('  Total CNS clearance (r_bc + r_bp):    %.4f h⁻¹\n', r_bc + r_bp);
fprintf('  Brain to CSF contribution:             %.2f%%\n', (r_bc / (r_bc + r_bp + r_cp)) * 100);
fprintf('  BBB contribution:                      %.2f%%\n\n', (r_bp / (r_bc + r_bp)) * 100);
fprintf('  Blood to CSF contribution:                      %.2f%%\n\n', (r_cp / (r_bc + r_bp)) * 100);

fprintf('SLEEP PERIOD:\n');
fprintf('  Brain to CSF efflux (r_bc × σ_bc):    %.4f h⁻¹\n', r_bc * sigma_bc);
fprintf('  Blood-Brain Barrier (r_bp × σ_bp):    %.4f h⁻¹\n', r_bp * sigma_bp);
fprintf('  Blood-CSF Barrier (r_cp × σ_cp):    %.4f h⁻¹\n', r_cp * sigma_cp);
fprintf('  Total CNS clearance:                   %.4f h⁻¹\n', r_bc * sigma_bc + r_bp * sigma_bp);
fprintf('  Brain to CSF contribution:             %.2f%%\n', (r_bc * sigma_bc / (r_bc * sigma_bc + r_bp * sigma_bp + r_cp * sigma_cp)) * 100);
fprintf('  BBB contribution:                      %.2f%%\n\n', (r_bp * sigma_bp / (r_bc * sigma_bc + r_bp * sigma_bp)) * 100);
fprintf('  Blood to CSF contribution:                      %.2f%%\n\n', (r_cp * sigma_cp / (r_bc * sigma_bc + r_bp * sigma_bp)) * 100);

% Plot
figure;
subplot(3, 1, 1);
hold on;
plot(t_100days_gwrmse1, sol_100days_gwrmse1(:, 1), 'k', 'LineWidth', 5.0, 'DisplayName', 'Brain Compartment1');
legend('show');
ylabel('Amyloid Concentration', 'FontSize', 12, 'FontWeight','bold');
xlim([2336, 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2360, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2376, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2384, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
hold off;

subplot(3, 1, 2);
hold on;
plot(t_100days_gwrmse1, sol_100days_gwrmse1(:, 2), 'k', 'LineWidth', 5.0, 'DisplayName', 'CSF Compartment1');
legend('show');
ylabel('Amyloid Concentration', 'FontSize', 12, 'FontWeight','bold');
xlim([2336, 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2360, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2376, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2384, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
hold off;

subplot(3, 1, 3);
hold on;
plot(t_100days_gwrmse1, sol_100days_gwrmse1(:, 3), 'k', 'LineWidth', 5.0, 'DisplayName', 'Plasma Compartment1');
legend('show');
xlabel('Time (hr)', 'FontSize', 12, 'FontWeight','bold');
ylabel('Amyloid Concentration', 'FontSize', 12, 'FontWeight','bold');
xlim([2336, 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2360, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2376, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
xline(2384, 'k--', 'LineWidth', 3, 'HandleVisibility', 'off', 'Alpha', 0.5);
hold off;
