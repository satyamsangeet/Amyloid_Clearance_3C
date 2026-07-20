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

% Baseline model: full 8hr sleep
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

% Partial sleep: 4h sleep only, then wake
function dydt_n = model2(t, y)
    r_bc = 0.038; r_bp = 0.014; r_cp = 0.00537;
    sigma_bc = 1.131; sigma_bp = 1.768; sigma_cp = 6.100;
    sigma_p = 4.253; A = 16.203; sigma_A = 0.772; r_p = 0.427;

    if (t >= 2352 && t < 2356)
        sw_cycle = 0;  % 4h sleep
    elseif (t >= 2356 && t < 2360)
        sw_cycle = 1;
    else
        sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    end

    dydt_n = zeros(3, 1);
    dydt_n(1) = A*sw_cycle + sigma_A*A*(1-sw_cycle) - (r_bp*sw_cycle + sigma_bp*r_bp*(1-sw_cycle) + r_bc*sw_cycle + sigma_bc*r_bc*(1-sw_cycle))*y(1);
    dydt_n(2) = (r_bc*sw_cycle + sigma_bc*r_bc*(1-sw_cycle))*y(1) - (r_cp*sw_cycle + sigma_cp*r_cp*(1-sw_cycle))*y(2);
    dydt_n(3) = (r_bp*sw_cycle + sigma_bp*r_bp*(1-sw_cycle))*y(1) + (r_cp*sw_cycle + sigma_cp*r_cp*(1-sw_cycle))*y(2) - (r_p*sw_cycle + sigma_p*r_p*(1-sw_cycle))*y(3);
end

% Fragmented sleep: 3x2h sleep bouts with 1h wake between = 6h total
function dydt_n = model3(t, y)
    r_bc = 0.038; r_bp = 0.014; r_cp = 0.00537;
    sigma_bc = 1.131; sigma_bp = 1.768; sigma_cp = 6.100;
    sigma_p = 4.253; A = 16.203; sigma_A = 0.772; r_p = 0.427;

    if (t >= 2352 && t < 2354) || (t >= 2355 && t < 2357) || (t >= 2358 && t < 2360)
        sw_cycle = 0;  % sleep (3 x 2h = 6h total)
    elseif (t >= 2354 && t < 2355) || (t >= 2357 && t < 2358)
        sw_cycle = 1;  % awakenings within sleep window (2 x 1h)
    else
        sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    end

    dydt_n = zeros(3, 1);
    dydt_n(1) = A*sw_cycle + sigma_A*A*(1-sw_cycle) - (r_bp*sw_cycle + sigma_bp*r_bp*(1-sw_cycle) + r_bc*sw_cycle + sigma_bc*r_bc*(1-sw_cycle))*y(1);
    dydt_n(2) = (r_bc*sw_cycle + sigma_bc*r_bc*(1-sw_cycle))*y(1) - (r_cp*sw_cycle + sigma_cp*r_cp*(1-sw_cycle))*y(2);
    dydt_n(3) = (r_bp*sw_cycle + sigma_bp*r_bp*(1-sw_cycle))*y(1) + (r_cp*sw_cycle + sigma_cp*r_cp*(1-sw_cycle))*y(2) - (r_p*sw_cycle + sigma_p*r_p*(1-sw_cycle))*y(3);
end

% Defining Euler-Maruyama Method
function [t, w] = euler(F, endpoints, initial_conditions, ts)
    total_time = endpoints(2) - endpoints(1);
    num_steps = floor(total_time / ts);
    
    t = linspace(endpoints(1), endpoints(2), num_steps + 1)';
    w = zeros(length(t), length(initial_conditions));
    w(1,:) = initial_conditions;

    for k = 1:num_steps
        w(k+1,:) = w(k,:) + F(t(k), w(k,:))' * ts;
    end
end

% Run simulation
[t_100days_gwrmse1, sol_100days_gwrmse1] = euler(@(t,y) model1(t,y), [0, 24*110], [0,600,20], 0.01);
[t_100days_gwrmse2, sol_100days_gwrmse2] = euler(@(t,y) model2(t,y), [0, 24*110], [0,600,20], 0.01);
[t_100days_gwrmse3, sol_100days_gwrmse3] = euler(@(t,y) model3(t,y), [0, 24*110], [0,600,20], 0.01);

x_start = 2312; x_end = 2456;
red_start = 2352; red_end = 2360;
green_start = 2360; green_end = 2456;

figure;
for sp = 1:3
    subplot(3, 1, sp);
    hold on;
    patch([red_start red_end red_end red_start], ...
          [-1e6 -1e6 1e6 1e6], [1 0.8 0.8], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
    patch([green_start green_end green_end green_start], ...
          [-1e6 -1e6 1e6 1e6], [0.8 1 0.8], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

    comp = sp;
    plot(t_100days_gwrmse1, sol_100days_gwrmse1(:, comp), 'k', 'LineWidth', 5.0, 'DisplayName', ['Normal - Comp ' num2str(comp)]);
    plot(t_100days_gwrmse2, sol_100days_gwrmse2(:, comp), 'r', 'LineWidth', 5.0, 'DisplayName', ['Partial SD - Comp ' num2str(comp)]);
    plot(t_100days_gwrmse3, sol_100days_gwrmse3(:, comp), 'b', 'LineWidth', 5.0, 'DisplayName', ['Fragmented SD - Comp ' num2str(comp)]);
    legend('show');
    ylabel('Amyloid Concentration', 'FontSize', 12, 'FontWeight', 'bold');

    xlim([x_start, x_end]);
    switch sp
        case 1
            ylim([220 300]);
        case 2
            ylim([650 900]);
        case 3
            ylim([15 20]);
    end

    tick_abs = x_start:4:x_end;
    tick_rel = tick_abs - x_start;
    xticks(tick_abs);
    xticklabels(arrayfun(@num2str, tick_rel, 'UniformOutput', false));
    xline(red_start,   'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);
    xline(red_end,     'k--', 'LineWidth', 2, 'HandleVisibility', 'off', 'Alpha', 0.5);

    if sp == 3
        xlabel('Time (hr)', 'FontSize', 12, 'FontWeight', 'bold');
    end
    hold off;
end
