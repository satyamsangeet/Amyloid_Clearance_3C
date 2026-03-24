% default graphuics
set(groot, ...
    'DefaultFigureColor', 'w', ...
    'DefaultAxesLineWidth', 1.2, ...
    'DefaultAxesXColor', 'k', ...
    'DefaultAxesYColor', 'k', ...
    'DefaultAxesFontUnits', 'points', ...
    'DefaultAxesFontSize', 10, ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultLineLineWidth', 1.5, ...
    'DefaultTextFontUnits', 'Points', ...
    'DefaultTextFontSize', 10, ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesBox', 'off', ...
    'DefaultAxesTickLength', [0.02 0.025]);

set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

function dydt_n = model_param_variation(t, y, params)
    A = params.A;
    sigma_A = params.sigma_A;
    r_bc = params.r_bc;
    sigma_bc = params.sigma_bc;
    r_cp = params.r_cp;
    sigma_cp = params.sigma_cp;
    r_bp = params.r_bp;
    sigma_bp = params.sigma_bp;
    r_p = params.r_p;
    sigma_p = params.sigma_p;
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);

    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function [t, w] = euler(F, endpoints, initial_conditions, ts, params)
    if length(endpoints) == 2
        h = ts;
        total_time = endpoints(2) - endpoints(1);
        num_steps = floor(total_time / h);
        t = linspace(endpoints(1), endpoints(2), num_steps + 1);
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints;
    end
    w = zeros(num_steps+1, length(initial_conditions));
    w(1,:) = initial_conditions;
    for k = 1:num_steps
        w(k+1,:) = w(k,:) + F(t(k), w(k,:), params)' * h;
    end
    t = t(:);
end

% Default parameters
params = struct();
params.A = 16.203;
params.sigma_A = 0.772;
params.r_bc = 0.038;
params.sigma_bc = 1.131;
params.r_cp = 0.00537; 
params.sigma_cp = 6.100;
params.r_bp = 0.014;
params.sigma_bp = 1.768;
params.r_p = 0.427;
params.sigma_p = 4.253;

% Generate variation
num_variations = 20;
default_value = 4.253; 

% array
sigma_p_values = linspace(0, 10, num_variations);
[~, default_index] = min(abs(sigma_p_values - default_value));
sigma_p_values(default_index) = default_value;

jet_colors = jet(num_variations);
colors = jet_colors;
colors(default_index,:) = [0 0 0];

figure('Position', [100 100 1200 800], ...
    'Renderer', 'painters', ...
    'Units', 'inches', ...
    'PaperPositionMode', 'auto');

titles = {'Brain Concentration', 'CSF Concentration', 'Plasma Concentration'};
ylabels = {'Brain Concentration', 'CSF Concentration', 'Plasma Concentration'};

for i = 1:3
    subplot(3, 1, i);
    hold on;
    
    for j = 1:length(sigma_p_values)
        params.sigma_p = sigma_p_values(j);
        [t, sol] = euler(@(t,y,p) model_param_variation(t,y,p), [0, 24*100], [1,600,15.5], 0.01, params);
        
        if j == default_index  % Default value
            plot(t(t >= 2336 & t <= 2386), sol(t >= 2336 & t <= 2386, i), ...
                 'Color', colors(j,:), 'LineWidth', 2.0);
        else
            plot(t(t >= 2336 & t <= 2386), sol(t >= 2336 & t <= 2386, i), ...
                 'Color', colors(j,:), 'LineWidth', 1.5);
        end
    end
    
    xline(2352, 'k--', 'LineWidth', 2.0, 'Alpha', 0.8);
    xline(2360, 'k--', 'LineWidth', 2.0, 'Alpha', 0.8);
    xline(2376, 'k--', 'LineWidth', 2.0, 'Alpha', 0.8);

    xlabel('Time (hr)', 'FontWeight', 'bold');
    ylabel(ylabels{i}, 'FontWeight', 'bold');
    title(titles{i}, 'FontWeight', 'bold', 'FontSize', 12);
    xlim([2336, 2384]);
    xticks(2336:2:2384);
    xticklabels(0:2:48);
    grid on;
    
    if i == 1
        legend_entries = cell(1, length(sigma_p_values));
        for j = 1:length(sigma_p_values)
            if j == default_index
                legend_entries{j} = sprintf('sigma_{p} = %.4f (default)', sigma_p_values(j));
            else
                legend_entries{j} = sprintf('sigma_{p} = %.4f', sigma_p_values(j));
            end
        end
        legend(legend_entries, 'Location', 'eastoutside', 'FontSize', 8);
    end
    
    hold off;
end

set(gcf, 'Units', 'inches');
set(gcf, 'PaperPositionMode', 'auto');
han = axes(gcf, 'visible', 'off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
sgtitle('Parameter Variation Analysis for sigma_{p}', 'FontSize', 14, 'FontWeight', 'bold');
%print(gcf, '/home/satyam/Documents/combined_fitting_data/final_fitting_wrmse_rmse_global_local/images/parameter_variation_analysis_sigma_cp.svg', '-dsvg', '-painters', '-r300');
