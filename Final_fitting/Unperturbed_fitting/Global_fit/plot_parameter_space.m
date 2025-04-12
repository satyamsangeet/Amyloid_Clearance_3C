function plot_parameter_space(solutions, num_starts)
    % Extracting params and obj values
    params_matrix = zeros(num_starts, 10);
    fvals = zeros(num_starts, 1);
    for i = 1:num_starts
        params_matrix(i,:) = solutions(i).params;
        fvals(i) = solutions(i).Fval;
    end

    param_names = {'r_bc', 'r_bp', 'r_cp', 'sigma_bc', 'sigma_bp', 'sigma_cp', 'sigma_p', 'A', 'sigma_A', 'r_p'};
    
    % colormap
    colors = jet(num_starts);
    
    figure('Position', [100, 100, 1200, 800]);

    for i = 1:10
        subplot(3, 4, i);
        hold on;
        
        for j = 1:num_starts
            scatter(params_matrix(j,i), fvals(j), 50, colors(j,:), 'filled');
        end
        
        xlabel(param_names{i});
        ylabel('Objective Value');
        title(sprintf('Effect of %s on Objective', param_names{i}));
        grid on;
        hold off;
    end

    subplot(3, 4, 11:12);
    hold on;
    
    h = zeros(num_starts, 1);
    for j = 1:num_starts
        h(j) = scatter(NaN, NaN, 50, colors(j,:), 'filled');
    end
    
    legend(h, arrayfun(@(x) sprintf('Run %d', x), 1:num_starts, 'UniformOutput', false), 'Location', 'bestoutside', 'NumColumns', ceil(num_starts/10));
    axis off;
    hold off;
    sgtitle('Parameter Space Exploration');
    
    % Optional
    %saveas(gcf, 'global_all_params/parameter_space.png');
    %print('global_all_params/parameter_space_highres', '-dpng', '-r300');
    
    % CSV save
    for param_idx = 1:10
        param_name = param_names{param_idx};
        csv_filename = sprintf('global_all_params/parameter_space_%s.csv', param_name);
        
        % col1=run no. | col2=param value | col3=obj value
        csv_data = zeros(num_starts, 3);
        for run_idx = 1:num_starts
            csv_data(run_idx, 1) = run_idx;
            csv_data(run_idx, 2) = params_matrix(run_idx, param_idx);
            csv_data(run_idx, 3) = fvals(run_idx);
        end
        
        writematrix(csv_data, csv_filename);
        fprintf('Saved parameter space data for %s to %s\n', param_name, csv_filename);
    end
end
