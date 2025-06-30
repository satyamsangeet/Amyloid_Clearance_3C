function total_error = objective_function(params, exp_csf1, exp_csf2, exp_csf3, exp_plasma1)
    global optimization_history;
    optimization_history = [optimization_history; params];

    % Sims
    [t, sol] = euler_solver(@(t,y) model(t,y,params), [0, 24*100], [0,600,15.5], 0.01);
    
    % Extracting 36 hours
    csf_data = sol(233600:237400,2);
    plasma_data = sol(233600:237400,3);
    
    % Selecting data points
    % Imp - should correspond to the experimnetal data points
    time_indices = 1:200:3801;
    csf_model = csf_data(time_indices);
    plasma_model = plasma_data(time_indices);
    
    % Selecting specific indices
    idx1 = [1:8, 14:20];
    idx2 = [1:9, 14:20];
    idx3 = [1:7, 14:20];
    csf_model1 = csf_model(idx1); csf_model2 = csf_model(idx2);
    csf_model3 = csf_model(idx3); plasma_model1 = plasma_model(idx3);
    
    % NRMSE
    total_error = calculate_nrmse(csf_model1, csf_model2, csf_model3, plasma_model1, exp_csf1, exp_csf2, exp_csf3, exp_plasma1);
    fprintf('Parameters: r_bc=%.4f, r_bp=%.4f, r_cp=%.4f, sigma_bc=%.4f, sigma_bp=%.4f, sigma_cp=%.4f, sigma_p=%.4f, A=%.4f, sigma_A=%.4f, r_p=%.4f\n', params);
    fprintf('Total NRMSE: %.4f\n', total_error);
end

function nrmse = calculate_nrmse(model1, model2, model3, model4, exp1, exp2, exp3, exp4)
    % Ensure column vectors
    model1 = model1(:);
    model2 = model2(:);
    model3 = model3(:);
    model4 = model4(:);
    exp1 = exp1(:);
    exp2 = exp2(:);
    exp3 = exp3(:);
    exp4 = exp4(:);
    
    % Calculate RMSE
    rmse1 = sqrt(mean((model1 - exp1).^2));
    rmse2 = sqrt(mean((model2 - exp2).^2));
    rmse3 = sqrt(mean((model3 - exp3).^2));
    rmse4 = sqrt(mean((model4 - exp4).^2));
    
    % Normalise
    nrmse1 = rmse1/abs(max(exp1) - min(exp1));
    nrmse2 = rmse2/abs(max(exp2) - min(exp2));
    nrmse3 = rmse3/abs(max(exp3) - min(exp3));
    nrmse4 = rmse4/abs(max(exp4) - min(exp4));
    
    nrmse = mean([nrmse1, nrmse2, nrmse3, nrmse4]);
end
