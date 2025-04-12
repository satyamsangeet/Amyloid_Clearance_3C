function [t, w] = euler_solver(F, endpoints, initial_conditions, ts)
    if length(endpoints) == 2
        h = ts;
        total_time = endpoints(2) - endpoints(1);
        num_steps = floor(total_time / h);
        t = linspace(endpoints(1), endpoints(2), num_steps + 1)';
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints(:);
    end
    
    w = zeros(num_steps+1, length(initial_conditions));
    w(1,:) = initial_conditions;
    
    for k = 1:num_steps
        w(k+1,:) = w(k,:) + F(t(k), w(k,:))' * h;
    end
end
