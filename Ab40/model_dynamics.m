% general graphics, this will apply to any figure you open
% (groot is the default figure object).
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
 
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

% Updated model function to optimize only sigma_bp (a), sigma_cp (b), and rbc (a12_wake)
function dydt_n = model(t, y)
    r_bc = 0.038;
    r_bp = 0.014;
    r_cp = 0.00537;
    sigma_bc = 1.131;
    sigma_bp = 1.768;
    sigma_cp = 6.100;
    sigma_p = 4.253;
    A = 16.203*8.81;
    sigma_A = 0.772;
    r_p = 0.427/1.44;

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

% Main script
% Run simulation
[t_100days_gwrmse1, sol_100days_gwrmse1] = euler(@(t,y) model(t,y), [0, 24*100], [800,600,20], 0.01);

% Parameters for clearance calculations
r_bc = 0.044;
r_bp = 0.020;
r_cp = 0.00593;
sigma_bc = 1.104;
sigma_bp = 2.412;
sigma_cp = 6.225;

% Brain volume (typical human brain ISF volume in mL)
V_brain = 280; % mL (interstitial fluid volume, can adjust based on your model)

% Extract steady-state concentrations from t=2336 to t=2384
analysis_idx = find(t_100days_gwrmse1 >= 2336 & t_100days_gwrmse1 <= 2384);
t_analysis = t_100days_gwrmse1(analysis_idx);
sol_analysis = sol_100days_gwrmse1(analysis_idx, :);

% Identify wake (8-24) and sleep (0-8) periods
time_of_day = mod(t_analysis, 24);
wake_idx = (time_of_day >= 8) & (time_of_day < 24);
sleep_idx = (time_of_day >= 0) & (time_of_day < 8);

% Calculate mean concentrations during wake and sleep
wake_brain_conc = mean(sol_analysis(wake_idx, 1));
sleep_brain_conc = mean(sol_analysis(sleep_idx, 1));
wake_csf_conc = mean(sol_analysis(wake_idx, 2));
sleep_csf_conc = mean(sol_analysis(sleep_idx, 2));

fprintf('\n=== CNS CLEARANCE ANALYSIS ===\n\n');

fprintf('STEADY-STATE BRAIN CONCENTRATIONS (t=2336 to t=2384):\n');
fprintf('  During wake:  %.2f pg/mL\n', wake_brain_conc);
fprintf('  During sleep: %.2f pg/mL\n\n', sleep_brain_conc);

fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('FRACTIONAL CLEARANCE RATE CONSTANTS (h⁻¹):\n');
fprintf('─────────────────────────────────────────────────────────────\n\n');

fprintf('WAKE PERIOD (8:00 - 24:00):\n');
fprintf('  Brain to CSF efflux (r_bc):           %.4f h⁻¹\n', r_bc);
fprintf('  Blood-Brain Barrier (r_bp):           %.4f h⁻¹\n', r_bp);
fprintf('  Blood-CSF Barrier (r_cp):           %.4f h⁻¹\n', r_cp);
fprintf('  Total CNS clearance (r_bc + r_bp):    %.4f h⁻¹\n', r_bc + r_bp + r_cp);
fprintf('  Brain to CSF contribution:             %.2f%%\n', (r_bc / (r_bc + r_bp + r_cp)) * 100);
fprintf('  BBB contribution:                      %.2f%%\n\n', (r_bp / (r_bc + r_bp + r_cp)) * 100);
fprintf('  Blood to CSF contribution:                      %.2f%%\n\n', (r_cp / (r_bc + r_bp + r_cp)) * 100);

fprintf('SLEEP PERIOD (0:00 - 8:00):\n');
fprintf('  Brain to CSF efflux (r_bc × σ_bc):    %.4f h⁻¹\n', r_bc * sigma_bc);
fprintf('  Blood-Brain Barrier (r_bp × σ_bp):    %.4f h⁻¹\n', r_bp * sigma_bp);
fprintf('  Blood-CSF Barrier (r_cp × σ_cp):    %.4f h⁻¹\n', r_cp * sigma_cp);
fprintf('  Total CNS clearance:                   %.4f h⁻¹\n', r_bc * sigma_bc + r_bp * sigma_bp + r_cp * sigma_cp);
fprintf('  Brain to CSF contribution:             %.2f%%\n', (r_bc * sigma_bc / (r_bc * sigma_bc + r_bp * sigma_bp + r_cp * sigma_cp)) * 100);
fprintf('  BBB contribution:                      %.2f%%\n\n', (r_bp * sigma_bp / (r_bc * sigma_bc + r_bp * sigma_bp + r_cp * sigma_cp)) * 100);
fprintf('  Blood to CSF contribution:                      %.2f%%\n\n', (r_cp * sigma_cp / (r_bc * sigma_bc + r_bp * sigma_bp + r_cp * sigma_cp)) * 100);

fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('ACTUAL FLUX CONTRIBUTIONS (rate × concentration):\n');
fprintf('─────────────────────────────────────────────────────────────\n\n');

% Calculate fluxes during wake (in model units: pg/mL × h⁻¹)
flux_bc_wake = r_bc * wake_brain_conc;
flux_bp_wake = r_bp * wake_brain_conc;
flux_cp_wake = r_cp * wake_brain_conc;
flux_total_wake = flux_bc_wake + flux_bp_wake + flux_cp_wake;

fprintf('WAKE PERIOD:\n');
fprintf('  Brain to CSF flux (r_bc × [Aβ]_brain):      %.2f pg/mL/h\n', flux_bc_wake);
fprintf('  BBB flux (r_bp × [Aβ]_brain):               %.2f pg/mL/h\n', flux_bp_wake);
fprintf('  CSF to Blood flux (r_cp × [Aβ]_CSF):               %.2f pg/mL/h\n', flux_cp_wake);
fprintf('  Total CNS clearance flux:                   %.2f pg/mL/h\n', flux_total_wake);
fprintf('  Brain to CSF flux contribution:             %.2f%%\n', (flux_bc_wake / flux_total_wake) * 100);
fprintf('  BBB flux contribution:                      %.2f%%\n\n', (flux_bp_wake / flux_total_wake) * 100);
fprintf('  Blood-CSF flux contribution:                      %.2f%%\n\n', (flux_cp_wake / flux_total_wake) * 100);

% Calculate fluxes during sleep
flux_bc_sleep = r_bc * sigma_bc * sleep_brain_conc;
flux_bp_sleep = r_bp * sigma_bp * sleep_brain_conc;
flux_cp_sleep = r_cp * sigma_cp * sleep_csf_conc;
flux_total_sleep = flux_bc_sleep + flux_bp_sleep + flux_cp_sleep;

fprintf('SLEEP PERIOD:\n');
fprintf('  Brain to CSF flux (r_bc × σ_bc × [Aβ]_brain):   %.2f pg/mL/h\n', flux_bc_sleep);
fprintf('  BBB flux (r_bp × σ_bp × [Aβ]_brain):            %.2f pg/mL/h\n', flux_bp_sleep);
fprintf('  Blood CSF flux (r_cp × σ_cp × [Aβ]_CSF):            %.2f pg/mL/h\n', flux_cp_sleep);
fprintf('  Total CNS clearance flux:                       %.2f pg/mL/h\n', flux_total_sleep);
fprintf('  Brain to CSF flux contribution:                 %.2f%%\n', (flux_bc_sleep / flux_total_sleep) * 100);
fprintf('  BBB flux contribution:                          %.2f%%\n\n', (flux_bp_sleep / flux_total_sleep) * 100);
fprintf('  Blood CSF flux contribution:                          %.2f%%\n\n', (flux_cp_sleep / flux_total_sleep) * 100);

fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('CONVERSION TO ng/min FOR COMPARISON WITH ROBERTS ET AL. (2014):\n');
fprintf('─────────────────────────────────────────────────────────────\n\n');

fprintf('Using brain ISF volume V_brain = %.0f mL\n\n', V_brain);

% Convert flux from pg/mL/h to ng/min
% flux (pg/mL/h) × V_brain (mL) × (1 ng / 1000 pg) × (1 h / 60 min) = ng/min
conversion_factor = V_brain * (1/1000) * (1/60);

flux_bc_wake_ngmin = flux_bc_wake * conversion_factor;
flux_bp_wake_ngmin = flux_bp_wake * conversion_factor;
flux_cp_wake_ngmin = flux_cp_wake * conversion_factor;
flux_total_wake_ngmin = flux_total_wake * conversion_factor;

flux_bc_sleep_ngmin = flux_bc_sleep * conversion_factor;
flux_bp_sleep_ngmin = flux_bp_sleep * conversion_factor;
flux_cp_sleep_ngmin = flux_cp_sleep * conversion_factor;
flux_total_sleep_ngmin = flux_total_sleep * conversion_factor;

fprintf('WAKE PERIOD (in ng/min):\n');
fprintf('  Brain to CSF flux:                      %.2f ng/min\n', flux_bc_wake_ngmin);
fprintf('  BBB flux:                               %.2f ng/min\n', flux_bp_wake_ngmin);
fprintf('  Blood CSF flux:                               %.2f ng/min\n', flux_cp_wake_ngmin);
fprintf('  Total CNS clearance flux:               %.2f ng/min\n', flux_total_wake_ngmin);
fprintf('  Brain to CSF contribution:              %.2f%%\n', (flux_bc_wake_ngmin / flux_total_wake_ngmin) * 100);
fprintf('  BBB contribution:                       %.2f%%\n\n', (flux_bp_wake_ngmin / flux_total_wake_ngmin) * 100);
fprintf('  Blood CSF contribution:                       %.2f%%\n\n', (flux_cp_wake_ngmin / flux_total_wake_ngmin) * 100);

fprintf('SLEEP PERIOD (in ng/min):\n');
fprintf('  Brain to CSF flux:                      %.2f ng/min\n', flux_bc_sleep_ngmin);
fprintf('  BBB flux:                               %.2f ng/min\n', flux_bp_sleep_ngmin);
fprintf('  Blood CSF flux:                               %.2f ng/min\n', flux_cp_sleep_ngmin);
fprintf('  Total CNS clearance flux:               %.2f ng/min\n', flux_total_sleep_ngmin);
fprintf('  Brain to CSF contribution:              %.2f%%\n', (flux_bc_sleep_ngmin / flux_total_sleep_ngmin) * 100);
fprintf('  BBB contribution:                       %.2f%%\n\n', (flux_bp_sleep_ngmin / flux_total_sleep_ngmin) * 100);
fprintf('  Blood CSF contribution:                       %.2f%%\n\n', (flux_cp_sleep_ngmin / flux_total_sleep_ngmin) * 100);

fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('COMPARISON WITH ROBERTS ET AL. (2014) EXPERIMENTAL DATA:\n');
fprintf('─────────────────────────────────────────────────────────────\n\n');

fprintf('Roberts et al. (2014) reported (presumably during wake):\n');
fprintf('  Total CNS clearance:    9.7 ng/min\n');
fprintf('  CSF pathway:            5.0 ng/min (51.5%%)\n');
fprintf('  BBB pathway:            4.7 ng/min (48.5%%)\n\n');

fprintf('Your model during WAKE:\n');
fprintf('  Total CNS clearance:    %.2f ng/min\n', flux_total_wake_ngmin);
fprintf('  CSF pathway:            %.2f ng/min (%.1f%%)\n', flux_bc_wake_ngmin, (flux_bc_wake_ngmin / flux_total_wake_ngmin) * 100);
fprintf('  BBB pathway:            %.2f ng/min (%.1f%%)\n\n', flux_bp_wake_ngmin, (flux_bp_wake_ngmin / flux_total_wake_ngmin) * 100);
fprintf('  Blood CSF pathway:            %.2f ng/min (%.1f%%)\n\n', flux_cp_wake_ngmin, (flux_cp_wake_ngmin / flux_total_wake_ngmin) * 100);

fprintf('Scaling factor needed to match Roberts et al. total flux:\n');
fprintf('  %.2fx (apply to concentration or volume if needed)\n\n', 9.7 / flux_total_wake_ngmin);

fprintf('─────────────────────────────────────────────────────────────\n');
fprintf('SLEEP vs WAKE COMPARISON:\n');
fprintf('─────────────────────────────────────────────────────────────\n\n');

fprintf('FLUX FOLD CHANGES (accounting for concentration):\n');
fprintf('  Brain to CSF flux increase:            %.2fx\n', flux_bc_sleep / flux_bc_wake);
fprintf('  BBB flux increase:                     %.2fx\n', flux_bp_sleep / flux_bp_wake);
fprintf('  Blood CSF flux increase:                     %.2fx\n', flux_cp_sleep / flux_cp_wake);
fprintf('  Total clearance flux increase:         %.2fx\n\n', flux_total_sleep / flux_total_wake);

% Plot comparison
figure;
subplot(3, 1, 1);
hold on;
plot(t_100days_gwrmse1, sol_100days_gwrmse1(:, 1), 'r', 'LineWidth', 3.0, 'DisplayName', 'Brain Compartment1');
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
plot(t_100days_gwrmse1, sol_100days_gwrmse1(:, 2), 'r', 'LineWidth', 3.0, 'DisplayName', 'CSF Compartment1');
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
plot(t_100days_gwrmse1, sol_100days_gwrmse1(:, 3), 'r', 'LineWidth', 3.0, 'DisplayName', 'Plasma Compartment1');
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

fprintf('═════════════════════════════════════════════════════════════\n');
fprintf('OSCILLATION AMPLITUDE CHANGES (WAKE → SLEEP):\n');
fprintf('═════════════════════════════════════════════════════════════\n\n');

% Define the two time points for comparison
t_wake_end = 2352;  % End of first wake period
t_sleep_end = 2360; % End of first sleep period

% Find indices closest to these times
[~, idx_wake] = min(abs(t_100days_gwrmse1 - t_wake_end));
[~, idx_sleep] = min(abs(t_100days_gwrmse1 - t_sleep_end));

% Extract concentrations at these time points
conc_wake = sol_100days_gwrmse1(idx_wake, :);
conc_sleep = sol_100days_gwrmse1(idx_sleep, :);

% Calculate percentage changes
pct_change_brain = ((conc_sleep(1) - conc_wake(1)) / conc_wake(1)) * 100;
pct_change_csf = ((conc_sleep(2) - conc_wake(2)) / conc_wake(2)) * 100;
pct_change_plasma = ((conc_sleep(3) - conc_wake(3)) / conc_wake(3)) * 100;

fprintf('Time comparison: t=%.0f (end of wake) → t=%.0f (end of sleep)\n\n', t_wake_end, t_sleep_end);

fprintf('BRAIN COMPARTMENT:\n');
fprintf('  Concentration at t=%.0f:  %.2f pg/mL\n', t_wake_end, conc_wake(1));
fprintf('  Concentration at t=%.0f:  %.2f pg/mL\n', t_sleep_end, conc_sleep(1));
fprintf('  Change:                    %+.2f%%\n\n', pct_change_brain);

fprintf('CSF COMPARTMENT:\n');
fprintf('  Concentration at t=%.0f:  %.2f pg/mL\n', t_wake_end, conc_wake(2));
fprintf('  Concentration at t=%.0f:  %.2f pg/mL\n', t_sleep_end, conc_sleep(2));
fprintf('  Change:                    %+.2f%%\n\n', pct_change_csf);

fprintf('PLASMA COMPARTMENT:\n');
fprintf('  Concentration at t=%.0f:  %.2f pg/mL\n', t_wake_end, conc_wake(3));
fprintf('  Concentration at t=%.0f:  %.2f pg/mL\n', t_sleep_end, conc_sleep(3));
fprintf('  Change:                    %+.2f%%\n\n', pct_change_plasma);

fprintf('SUMMARY:\n');
if pct_change_brain < 0
    fprintf('  Brain:  %.2f%% DECREASE during sleep\n', abs(pct_change_brain));
else
    fprintf('  Brain:  %.2f%% INCREASE during sleep\n', pct_change_brain);
end

if pct_change_csf < 0
    fprintf('  CSF:    %.2f%% DECREASE during sleep\n', abs(pct_change_csf));
else
    fprintf('  CSF:    %.2f%% INCREASE during sleep\n', pct_change_csf);
end

if pct_change_plasma < 0
    fprintf('  Plasma: %.2f%% DECREASE during sleep\n\n', abs(pct_change_plasma));
else
    fprintf('  Plasma: %.2f%% INCREASE during sleep\n\n', pct_change_plasma);
end