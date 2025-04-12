% Function to plot parameter space
function plot_parameter_space(solutions, num_starts)
    % Extract parameters and objective values
    params_matrix = zeros(num_starts, 10);
    fvals = zeros(num_starts, 1);
    for i = 1:num_starts
        params_matrix(i,:) = solutions(i).params;
        fvals(i) = solutions(i).Fval;
    end
    
    % Create parameter names for plotting
    param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
    
    % Get colormap with num_starts distinct colors
    colors = jet(num_starts);
    
    % Create figure with subplots for each parameter
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot parameters in subplots
    for i = 1:10
        subplot(3, 4, i);
        hold on;
        
        % Plot each run's data point with a distinct color
        for j = 1:num_starts
            scatter(params_matrix(j,i), fvals(j), 50, colors(j,:), 'filled');
        end
        
        xlabel(param_names{i});
        ylabel('Objective Value');
        title(sprintf('Effect of %s on Objective', param_names{i}));
        grid on;
        hold off;
    end
    
    % Create a separate subplot for the legend
    subplot(3, 4, 11:12);
    hold on;
    
    % Create dummy points for legend
    h = zeros(num_starts, 1);
    for j = 1:num_starts
        h(j) = scatter(NaN, NaN, 50, colors(j,:), 'filled');
    end
    
    % Add legend with run numbers
    legend(h, arrayfun(@(x) sprintf('Run %d', x), 1:num_starts, 'UniformOutput', false), 'Location', 'bestoutside', 'NumColumns', ceil(num_starts/10));
    axis off;
    hold off;
    
    % Add title to the figure
    sgtitle('Parameter Space Exploration');
    
    % Save the figure
    saveas(gcf, 'global_all_params/parameter_space.png');
    print('global_all_params/parameter_space_highres', '-dpng', '-r300');
    
    % NEW CODE: Save parameter space data to CSV files
    for param_idx = 1:10
        % Create and open CSV file for this parameter
        param_name = param_names{param_idx};
        csv_filename = sprintf('global_all_params/parameter_space_%s.csv', param_name);
        
        % Prepare data in the required format: [run_number, parameter_value, objective_value]
        csv_data = zeros(num_starts, 3);
        for run_idx = 1:num_starts
            csv_data(run_idx, 1) = run_idx;  % Run number
            csv_data(run_idx, 2) = params_matrix(run_idx, param_idx);  % Parameter value
            csv_data(run_idx, 3) = fvals(run_idx);  % Objective function value
        end
        
        % Save to CSV
        writematrix(csv_data, csv_filename);
        fprintf('Saved parameter space data for %s to %s\n', param_name, csv_filename);
    end
end
