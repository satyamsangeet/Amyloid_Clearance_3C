function [solutions, best_params, best_fval, all_trajectories] = run_optimization(problem, start_points, num_starts, data)
    global optimization_history;

    solutions = struct('params', cell(1, num_starts), 'Fval', cell(1, num_starts), ...
        'Exitflag', cell(1, num_starts));
    all_trajectories = cell(num_starts, 1);
    best_fval = Inf;
    best_params = [];
    
    if ~exist('global_all_params', 'dir')
        mkdir('global_all_params');
    end
    
    for i = 1:num_starts
        optimization_history = [];
        problem.x0 = start_points(i,:);
        
        [x, fval, exitflag] = fmincon(problem);
        solutions(i).params = x(:)';
        solutions(i).Fval = fval;
        solutions(i).Exitflag = exitflag;
        all_trajectories{i} = optimization_history;
        
        if fval < best_fval
            best_fval = fval;
            best_params = x;
        end
        
        save(sprintf('global_all_params/run_%d_results.mat', i), 'x', 'fval', 'exitflag', 'optimization_history');
    end
    
    save('global_all_params/final_results.mat', 'solutions', 'best_params', 'best_fval', 'all_trajectories');
end

function start_points = generate_start_points(initial_guess, lb, ub, num_starts)
    start_points = zeros(num_starts, length(initial_guess));
    start_points(1,:) = initial_guess;
    for i = 2:num_starts
        start_points(i,:) = lb + rand(1, length(lb)) .* (ub - lb);
    end
end
