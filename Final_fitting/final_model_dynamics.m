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
function dydt_n = model_global_wrmse(t, y)
    A = 12;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_cp = 0.0056;
    sigma_cp = 3.99;
    r_bp = (r_bc*(1-133*r_cp))/(133*r_cp);
    sigma_bp = 3.99;
    r_p = 0.28;
    sigma_p = 2.89;

    % Switch between sleep and wake states
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % Define the ODE system
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
[t_100days_gwrmse, sol_100days_gwrmse] = euler(@(t,y) model_global_wrmse(t,y), [0, 24*100], [0,600,15.5], 0.01);

% Plot comparison
figure;
subplot(3, 1, 1);
hold on;
plot(t_100days_gwrmse, sol_100days_gwrmse(:, 1), 'r', 'LineWidth', 3.0, 'DisplayName', 'Brain Compartment');
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
plot(t_100days_gwrmse, sol_100days_gwrmse(:, 2), 'r', 'LineWidth', 3.0, 'DisplayName', 'CSF Compartment');
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
plot(t_100days_gwrmse, sol_100days_gwrmse(:, 3), 'r', 'LineWidth', 3.0, 'DisplayName', 'Plasma Compartment');
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