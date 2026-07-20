clear; clc;
N      = 500;     % Base sample size. Total ODE solves = N*(2D+2)
tspan  = [0, 2400]; % 100 days
dt     = 0.1;       
y0     = [0, 600, 16];
params = {'A','sigma_A','r_bc','sigma_bc','r_bp','sigma_bp','r_cp','sigma_cp','r_p','sigma_p'};
D      = length(params);

param_bounds = [
     0.0,  111.0;   % A
     0.0,    0.99;  % sigma_A
     0.01,   0.25;  % r_bc
     1.0,    7.0;   % sigma_bc
     0.01,   0.25;  % r_bp
     1.0,    7.0;   % sigma_bp
     0.0,    0.01;  % r_cp
     1.0,    7.0;   % sigma_cp
     0.0,    0.6;   % r_p
     1.0,    7.0    % sigma_p
];

fprintf('Generating Saltelli sample matrices (N=%d, D=%d)...\n', N, D);
rng_sobol = sobolset(2*D, 'Skip', 1000, 'Leap', 100);
raw       = net(rng_sobol, N);   % N x 2D in [0,1]

raw_A = raw(:, 1:D);
raw_B = raw(:, D+1:2*D);

A_mat = zeros(N, D);
B_mat = zeros(N, D);
for k = 1:D
    lo = param_bounds(k, 1);
    hi = param_bounds(k, 2);
    A_mat(:, k) = lo + (hi - lo) * raw_A(:, k);
    B_mat(:, k) = lo + (hi - lo) * raw_B(:, k);
end

% Hybrid matrices
AB = cell(D, 1);
for i = 1:D
    M       = A_mat;       % copy of A
    M(:, i) = B_mat(:, i); % replace column i with B
    AB{i}   = M;
end

fprintf('Total ODE solves: %d\n\n', N * (2 + D));
fprintf('Evaluating model on matrix A...\n');
Y_A  = evaluate_model(A_mat, N, tspan, dt, y0);

fprintf('Evaluating model on matrix B...\n');
Y_B  = evaluate_model(B_mat, N, tspan, dt, y0);

Y_AB = zeros(N, D);
for i = 1:D
    fprintf('Evaluating model on AB_%d (%s)...\n', i, params{i});
    Y_AB(:, i) = evaluate_model(AB{i}, N, tspan, dt, y0);
end

valid = isfinite(Y_A) & isfinite(Y_B) & all(isfinite(Y_AB), 2);
fprintf('\nValid samples: %d / %d\n', sum(valid), N);

Y_A  = Y_A(valid);
Y_B  = Y_B(valid);
Y_AB = Y_AB(valid, :);
n    = sum(valid);

Y_all = [Y_A; Y_B; Y_AB(:)];
f0    = mean(Y_all);
V_Y   = var(Y_all);

fprintf('Output variance: %.4f\n', V_Y);

Si  = zeros(1, D);
STi = zeros(1, D);

for i = 1:D
    Si(i)  = (1/n) * sum(Y_B .* (Y_AB(:,i) - Y_A)) / V_Y;
    STi(i) = (1/(2*n)) * sum((Y_A - Y_AB(:,i)).^2) / V_Y;
end

Si  = max(0, min(1, Si));
STi = max(0, min(1, STi));

fprintf('%-12s  %12s  %12s  %12s\n', 'Parameter', 'First-order', 'Total-order', 'Interaction');
fprintf('%s\n', repmat('-', 1, 54));
for i = 1:D
    interaction = max(0, STi(i) - Si(i));
    fprintf('%-12s  %12.4f  %12.4f  %12.4f\n', params{i}, Si(i), STi(i), interaction);
end

fprintf('\nSum of first-order indices: %.4f  (should be <= 1)\n', sum(Si));
fprintf('Sum of total-order indices: %.4f\n\n', sum(STi));

% Ranking
fprintf('Parameter ranking by First-Order index:\n');
[sorted_Si, ord] = sort(Si, 'descend');
for i = 1:D
    fprintf('  %2d. %-12s  %.4f\n', i, params{ord(i)}, sorted_Si(i));
end
fprintf('\nParameter ranking by Total-Order index:\n');
[sorted_STi, ord2] = sort(STi, 'descend');
for i = 1:D
    fprintf('  %2d. %-12s  %.4f\n', i, params{ord2(i)}, sorted_STi(i));
end


results_table = table(params', Si', STi', max(0, STi'-Si'), 'VariableNames', {'Parameter','FirstOrder','TotalOrder','Interaction'});
writetable(results_table, 'sobol_results.csv');
fprintf('\nNumerical results saved to sobol_results.csv\n');

figure('Position', [100 100 900 500]);
colors = [0.18 0.44 0.71;   % First-order (blue)
          0.89 0.10 0.11];  % Total-order  (red)

subplot(1,2,1);
x = 1:D;
bar_data = [Si; STi]';
b = bar(x, bar_data, 'grouped');
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
set(gca, 'XTick', x, 'XTickLabel', params, 'XTickLabelRotation', 45, 'FontSize', 9, 'Box', 'off');
ylabel('Sensitivity Index', 'FontSize', 10);
title('Sobol Sensitivity Indices', 'FontSize', 11, 'FontWeight', 'bold');
legend({'First-order (S_i)', 'Total-order (S_{Ti})'}, 'Location', 'northwest', 'FontSize', 8);
grid on; grid minor;
ylim([0 1]);

subplot(1,2,2);
bar_stack = [Si; max(0, STi - Si)]';
b2 = bar(x, bar_stack, 'stacked');
b2(1).FaceColor = colors(1,:); b2(1).EdgeColor = 'none';
b2(2).FaceColor = [0.4 0.7 0.4]; b2(2).EdgeColor = 'none';
set(gca, 'XTick', x, 'XTickLabel', params, 'XTickLabelRotation', 45, 'FontSize', 9, 'Box', 'off');
ylabel('Sensitivity Index', 'FontSize', 10);
title('First-Order + Interaction Effects', 'FontSize', 11, 'FontWeight', 'bold');
legend({'First-order (S_i)', 'Interaction (S_{Ti}-S_i)'}, 'Location', 'northwest', 'FontSize', 8);
grid on; grid minor;
ylim([0, max(1, max(STi)*1.05)]);

sgtitle(sprintf('Saltelli-Sobol Analysis  (N=%d, %d ODE solves)', N, N*(2+D)), ...
        'FontSize', 12, 'FontWeight', 'bold');

saveas(gcf, 'sobol_sensitivity_correct.png');
fprintf('Plot saved to sobol_sensitivity_correct.png\n');

function Y = evaluate_model(samples, N, tspan, dt, y0)
    Y = NaN(N, 1);
    for i = 1:N
        A        = samples(i,1);
        sigma_A  = samples(i,2);
        r_bc     = samples(i,3);
        sigma_bc = samples(i,4);
        r_bp     = samples(i,5);
        sigma_bp = samples(i,6);
        r_cp     = samples(i,7);
        sigma_cp = samples(i,8);
        r_p      = samples(i,9);
        sigma_p  = samples(i,10);

        try
            [t, y] = euler( @(t,y) model3(t, y, A, sigma_A, r_bc, sigma_bc, r_cp, sigma_cp, r_bp, sigma_bp, r_p, sigma_p), tspan, y0, dt);

            final_cycle = t >= (tspan(2) - 24);
            if sum(final_cycle) > 0
                Y(i) = mean(y(final_cycle, 2));
            else
                Y(i) = y(end, 1);
            end

            if Y(i) < 0 || Y(i) > 1e8
                Y(i) = NaN;
            end
        catch
            Y(i) = NaN;
        end
    end
end
