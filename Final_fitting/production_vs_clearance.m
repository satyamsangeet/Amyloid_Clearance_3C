% Modified model function to test different hypotheses
function dydt_n = model_hypothesis(t, y, sigma_A, sigma_bc)
    % Base parameters
    A = 13;
    r_bc = 1.5;
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
    
    dydt_n(1) = (A * sw_cycle + sigma_A * A * (1 - sw_cycle)) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
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

% Run simulations for each hypothesis
[t, y1] = euler(@(t,y) model_hypothesis(t,y,0.2,1.0), [0,2400], [0,600,15.5], 0.01);  % Production
[~, y2] = euler(@(t,y) model_hypothesis(t,y,1.0,4.0), [0,2400], [0,600,15.5], 0.01);  % Clearance
[~, y3] = euler(@(t,y) model_hypothesis(t,y,0.2,4.0), [0,2400], [0,600,15.5], 0.01);  % Combined
[~, y4] = euler(@(t,y) model_hypothesis(t,y,0.8,2.5), [0,2400], [0,600,15.5], 0.01);  % Default

% Plot results
figure('Position', [100 100 800 600]);
subplot(3,1,1)
plot(t, y1(:,1), 'r-', 'LineWidth', 2, 'DisplayName', 'Production Only');
hold on;
plot(t, y2(:,1), 'b-', 'LineWidth', 2, 'DisplayName', 'Clearance Only');
plot(t, y3(:,1), 'g-', 'LineWidth', 2, 'DisplayName', 'Combined');
plot(t, y4(:,1), 'm-', 'LineWidth', 2, 'DisplayName', 'Default');
ylabel('Brain Amyloid Concentration', 'FontWeight', 'bold');
legend('show');
xlim([2336 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2360, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2376, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

subplot(3,1,2)
plot(t, y1(:,2), 'r-', 'LineWidth', 2, 'DisplayName', 'Production Only');
hold on;
plot(t, y2(:,2), 'b-', 'LineWidth', 2, 'DisplayName', 'Clearance Only');
plot(t, y3(:,2), 'g-', 'LineWidth', 2, 'DisplayName', 'Combined');
plot(t, y4(:,2), 'm-', 'LineWidth', 2, 'DisplayName', 'Default');
ylabel('CSF Amyloid Concentration', 'FontWeight', 'bold');
xlim([2336 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2360, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2376, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

subplot(3,1,3)
plot(t, y1(:,3), 'r-', 'LineWidth', 2, 'DisplayName', 'Production Only');
hold on;
plot(t, y2(:,3), 'b-', 'LineWidth', 2, 'DisplayName', 'Clearance Only');
plot(t, y3(:,3), 'g-', 'LineWidth', 2, 'DisplayName', 'Combined');
plot(t, y4(:,3), 'm-', 'LineWidth', 2, 'DisplayName', 'Default');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Plasma Amyloid Concentration', 'FontWeight', 'bold');
xlim([2336 2384]);
xticks(2336:2:2384);
xticklabels(0:2:48);
xline(2352, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2360, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xline(2376, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');