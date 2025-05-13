%H1
function dydt_n = model1(t, y)
    A = 12.330;
    sigma_A = 0.787;
    r_bp = 0.655;
    sigma_bp = 6.152;
    r_p = 0.337;
    sigma_p = 4.154;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        r_bc = 2.405;
        r_cp = 0.00458;
        sigma_bc = 1.193;
        sigma_cp = 5.539;
    else
        r_bc = 2.361;
        sigma_bc = 2.519;
        r_cp = 0.00458;
        sigma_cp = 5.540;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H2
function dydt_n = model2(t, y)
    A = 12.330;
    r_bp = 0.655;
    r_p = 0.337;
    r_bc = 2.361;
    r_cp = 0.00458;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        sigma_A = 0.791;
        sigma_bc = 2.479;
        sigma_bp = 5.890;
        sigma_cp = 4.864;
        sigma_p = 3.531;
    else
        sigma_A = 0.787;
        sigma_bc = 2.519;
        sigma_bp = 6.152;
        sigma_cp = 5.540;
        sigma_p = 4.154;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);

    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H3
function dydt_n = model3(t, y)
    A = 12.330;
    r_bp = 0.655;
    r_p = 0.337;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
    	r_bc = 2.421;
    	r_cp = 0.00433;
        sigma_A = 0.849;
        sigma_bc = 1.690;
        sigma_bp = 5.473;
        sigma_cp = 2.561;
        sigma_p = 2.366;
    else
    	r_bc = 2.361;
    	r_cp = 0.00458;
        sigma_A = 0.787;
        sigma_bc = 2.519;
        sigma_bp = 6.152;
        sigma_cp = 5.540;
        sigma_p = 4.154;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H4
function dydt_n = model4(t, y)
    A = 12.330;
    r_p = 0.337;
    sigma_A = 0.787;
    sigma_p = 4.154;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
    	r_bc = 4.086;
    	r_bp = 1.304;
    	r_cp = 0.00407;
        sigma_bc = 0.332;
        sigma_bp = 4.604;
        sigma_cp = 5.378;
    else
    	r_bc = 2.361;
    	r_cp = 0.00458;
        r_bp = 0.655;
        sigma_bc = 2.519;
        sigma_bp = 6.152;
        sigma_cp = 5.540;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H5
function dydt_n = model5(t, y)
    A = 12.330;
    sigma_A = 0.787;
    r_bp = 0.655;
    sigma_bp = 6.152;
    r_p = 0.337;
    sigma_p = 4.154;
    r_bc = 2.361;
    r_cp = 0.00458;

    % Switch between sleep and wake states
    if(t>=2352 && t<2372)
    	sigma_bc = 2.519;
        sigma_cp = 5.540;
    else
    	sigma_bc = 2.519;
        sigma_cp = 5.540;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);

    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%EM
function [t, w] = euler(F, endpoints, initial_conditions, ts)
    if length(endpoints) == 2
        h = ts; %delta_t (seconds)
        total_time = endpoints(2) - endpoints(1);
        num_steps = floor(total_time / h);
        t = linspace(endpoints(1), endpoints(2), num_steps + 1); %Creating time vector
    else
        h = endpoints(2) - endpoints(1);
        t = endpoints;
    end
    w = zeros(num_steps+1, length(initial_conditions));
    w(1,:) = initial_conditions;
    for k = 1:num_steps
        w(k+1,:) = w(k,:) + F(t(k), w(k,:))' * h;
    end
    t = t(:);
end

[t_100days_global1, sol_100days_global1] = euler(@(t,y) model1(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_global2, sol_100days_global2] = euler(@(t,y) model2(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_global3, sol_100days_global3] = euler(@(t,y) model3(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_global4, sol_100days_global4] = euler(@(t,y) model4(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_global5, sol_100days_global5] = euler(@(t,y) model5(t,y), [0, 24*100], [0,600,15.5], 0.01);

cdata_36hours1_global1 = sol_100days_global1(233600:237200, 2);
cdata_36hours1_global2 = sol_100days_global2(233600:237200, 2);
cdata_36hours1_global3 = sol_100days_global3(233600:237200, 2);
cdata_36hours1_global4 = sol_100days_global4(233600:237200, 2);
cdata_36hours1_global5 = sol_100days_global5(233600:237200, 2);

pdata_36hours1_global1 = sol_100days_global1(233600:237200, 3);
pdata_36hours1_global2 = sol_100days_global2(233600:237200, 3);
pdata_36hours1_global3 = sol_100days_global3(233600:237200, 3);
pdata_36hours1_global4 = sol_100days_global4(233600:237200, 3);
pdata_36hours1_global5 = sol_100days_global5(233600:237200, 3);

csf_data_file1 = 'data/liu2022_csf_concentration.csv';
plasma_data_file1 = 'data/liu2022_plasma_concentration.csv';
csf_data1 = readtable(csf_data_file1);
plasma_data1 = readtable(plasma_data_file1);

time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Concentration;
plasma_conc_exp1 = plasma_data1.Concentration;
csf_lsd1 = csf_data1.LSD;
plasma_lsd1 = plasma_data1.LSD;
csf_usd1 = csf_data1.USD;
plasma_usd1 = plasma_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;
plasma_std1 = (plasma_usd1 - plasma_lsd1)/2;

exp_csf1 = csf_conc_exp1(:);
exp_plasma1 = plasma_conc_exp1(:);

time_indices = 1:200:3601;

%NRMSE
function err = calculate_wrmse(model_data, exp_data)
    errors = (model_data - exp_data).^2;
    error = sqrt(sum(errors) / length(exp_data));
    err = error/abs(max(exp_data) - min(exp_data));
end

csf_36hr_model_global1 = cdata_36hours1_global1(time_indices);
csf_36hr_model_global2 = cdata_36hours1_global2(time_indices);
csf_36hr_model_global3 = cdata_36hours1_global3(time_indices);
csf_36hr_model_global4 = cdata_36hours1_global4(time_indices);
csf_36hr_model_global5 = cdata_36hours1_global5(time_indices);
plasma_36hr_model_global1 = pdata_36hours1_global1(time_indices);
plasma_36hr_model_global2 = pdata_36hours1_global2(time_indices);
plasma_36hr_model_global3 = pdata_36hours1_global3(time_indices);
plasma_36hr_model_global4 = pdata_36hours1_global4(time_indices);
plasma_36hr_model_global5 = pdata_36hours1_global5(time_indices);

wrmse1_global1 = calculate_wrmse(csf_36hr_model_global1, exp_csf1);
wrmse1_global2 = calculate_wrmse(csf_36hr_model_global2, exp_csf1);
wrmse1_global3 = calculate_wrmse(csf_36hr_model_global3, exp_csf1);
wrmse1_global4 = calculate_wrmse(csf_36hr_model_global4, exp_csf1);
wrmse1_global5 = calculate_wrmse(csf_36hr_model_global5, exp_csf1);
wrmse1_global6 = calculate_wrmse(plasma_36hr_model_global1, exp_plasma1);
wrmse1_global7 = calculate_wrmse(plasma_36hr_model_global2, exp_plasma1);
wrmse1_global8 = calculate_wrmse(plasma_36hr_model_global3, exp_plasma1);
wrmse1_global9 = calculate_wrmse(plasma_36hr_model_global4, exp_plasma1);
wrmse1_global10 = calculate_wrmse(plasma_36hr_model_global5, exp_plasma1);

new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);
colormap_jet = jet(5);

figure();
x1 = [time_exp1/10; flipud(time_exp1/10)];
inBetween1 = [csf_lsd1; flipud(csf_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('H1 (NRMSE: %.3f)', wrmse1_global1));
plot(x, cdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('H2 (NRMSE: %.3f)', wrmse1_global2));
plot(x, cdata_36hours1_global3, 'LineWidth', 2.0, 'Color', colormap_jet(3,:), 'DisplayName', sprintf('H3 (NRMSE: %.3f)', wrmse1_global3));
plot(x, cdata_36hours1_global4, 'LineWidth', 2.0, 'Color', colormap_jet(4,:), 'DisplayName', sprintf('H4 (NRMSE: %.3f)', wrmse1_global4));
plot(x, cdata_36hours1_global5, 'LineWidth', 2.0, 'Color', colormap_jet(5,:), 'DisplayName', sprintf('H5 (NRMSE: %.3f)', wrmse1_global5));
plot(time_exp1/10, csf_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Liu2022 - CSF');
model_points_global1 = interp1(x, cdata_36hours1_global1, time_exp1/10);
model_points_global2 = interp1(x, cdata_36hours1_global2, time_exp1/10);
errorbar(time_exp1/10, csf_conc_exp1, csf_std1, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1/10, csf_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1/10, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1/10, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2372]);
xticks(0:200:3600);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 3600]);
hold off;

figure();
x4 = [time_exp1/10; flipud(time_exp1/10)];
inBetween4 = [plasma_lsd1; flipud(plasma_usd1)];
fill(x4, inBetween4, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, pdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('H1 (NRMSE: %.3f)', wrmse1_global6));
plot(x, pdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('H2 (NRMSE: %.3f)', wrmse1_global7));
plot(x, pdata_36hours1_global3, 'LineWidth', 2.0, 'Color', colormap_jet(3,:), 'DisplayName', sprintf('H3 (NRMSE: %.3f)', wrmse1_global8));
plot(x, pdata_36hours1_global4, 'LineWidth', 2.0, 'Color', colormap_jet(4,:), 'DisplayName', sprintf('H4 (NRMSE: %.3f)', wrmse1_global9));
plot(x, pdata_36hours1_global5, 'LineWidth', 2.0, 'Color', colormap_jet(5,:), 'DisplayName', sprintf('H5 (NRMSE: %.3f)', wrmse1_global10));
plot(time_exp1/10, plasma_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global4 = interp1(x, pdata_36hours1_global1, time_exp1/10);
model_points_global5 = interp1(x, pdata_36hours1_global2, time_exp1/10);
errorbar(time_exp1/10, plasma_conc_exp1, plasma_std1, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1/10, plasma_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1/10, model_points_global4, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1/10, model_points_global5, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2372]);
xticks(0:200:3600);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Amyloid Concentration (pg/ml)', 'FontWeight', 'bold');
xlim([0, 3600]);
hold off;
