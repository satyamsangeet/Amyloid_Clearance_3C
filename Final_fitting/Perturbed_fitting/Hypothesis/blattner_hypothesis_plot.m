%H1
function dydt_n = model1(t, y)
    A = 13.791866842295299;
    sigma_A = 0.718493031801368;
    r_bp = 0.092859916742039;
    sigma_bp = 2.526995585282243;
    r_p = 0.285014821392297;
    sigma_p = 5.859219697739381;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        r_bc = 4.466;
        r_cp = 0.00596;
        sigma_bc = 0.073;
        sigma_cp = 0.920;
    else
        r_bc = 2.349030412757668;
        sigma_bc = 1.117763182229140;
        r_cp = 0.006556023403031;
        sigma_cp = 6.416831718871831;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H2
function dydt_n = model2(t, y)
    A = 13.791866842295299; 
    r_bp = 0.092859916742039;
    r_p = 0.285014821392297;
    r_bc = 2.349030412757668;
    r_cp = 0.006556023403031;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
        sigma_A = 0.899;
        sigma_bc = 0.223;
        sigma_bp = 0.100;
        sigma_cp = 1.557;
        sigma_p = 3.966;
    else
        sigma_A = 0.718493031801368;
        sigma_bc = 1.117763182229140;
        sigma_bp = 2.526995585282243;
        sigma_cp = 6.416831718871831;
        sigma_p = 5.859219697739381;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);

    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H3
function dydt_n = model3(t, y)
    A = 13.791866842295299; 
    r_bp = 0.092859916742039;
    r_p = 0.285014821392297;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
    	r_bc = 4.296;
    	r_cp = 0.0057;
        sigma_A = 0.865;
        sigma_bc = 0.123;
        sigma_bp = 1.189;
        sigma_cp = 1.921;
        sigma_p = 0.862;
    else
    	r_bc = 2.349030412757668;
    	r_cp = 0.006556023403031;
        sigma_A = 0.718493031801368;
        sigma_bc = 1.117763182229140;
        sigma_bp = 2.526995585282243;
        sigma_cp = 6.416831718871831;
        sigma_p = 5.859219697739381;
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H4
function dydt_n = model4(t, y)
    A = 13.791866842295299; 
    r_p = 0.285014821392297;
    sigma_A = 0.718493031801368;
    sigma_p = 5.859219697739381;

    % Switch between sleep and wake states
    if(t>=2336 && t<2372)
    	r_bc = 3.940;
    	r_bp = 1.681;
    	r_cp = 0.00047;
        sigma_bc = 0.733;
        sigma_bp = 2.844;
        sigma_cp = 4.441;
    else
    	r_bc = 2.349030412757668;
    	r_cp = 0.006556023403031;
        r_bp = 0.092859916742039;
        sigma_bc = 1.117763182229140;
        sigma_bp = 2.526995585282243;
        sigma_cp = 6.416831718871831;
        
    end

    sw_cycle = (mod(t,24) >=8 && mod(t,24) < 24);
    
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

%H5
function dydt_n = model5(t, y)
    A = 13.791866842295299;
    sigma_A = 0.718493031801368;
    r_bp = 0.092859916742039;
    sigma_bp = 2.526995585282243;
    r_p = 0.285014821392297;
    sigma_p = 5.859219697739381;
    r_cp = 0.006556023403031;
    r_bc = 2.349030412757668;

    % Switch between sleep and wake states
    if(t>=2352 && t<2372)
        sigma_bc = 0.1003;
        sigma_cp = 0.2574;
    else
        sigma_bc = 1.117763182229140;
        sigma_cp = 6.416831718871831;
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

csf_data_file1 = 'data/blattner2020_csf_concentration.csv';
csf_data1 = readtable(csf_data_file1);

time_exp1 = csf_data1.Time;
csf_conc_exp1 = csf_data1.Concentration;
csf_lsd1 = csf_data1.LSD;
csf_usd1 = csf_data1.USD;
csf_std1 = (csf_usd1 - csf_lsd1)/2;
exp_csf1 = csf_conc_exp1(:);

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

wrmse1_global1 = calculate_wrmse(csf_36hr_model_global1, exp_csf1);
wrmse1_global2 = calculate_wrmse(csf_36hr_model_global2, exp_csf1);
wrmse1_global3 = calculate_wrmse(csf_36hr_model_global3, exp_csf1);
wrmse1_global4 = calculate_wrmse(csf_36hr_model_global4, exp_csf1);
wrmse1_global5 = calculate_wrmse(csf_36hr_model_global5, exp_csf1);

new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);
colormap_jet = jet(5);

x = 1:3601;

figure();
x1 = [time_exp1; flipud(time_exp1)];
inBetween1 = [csf_lsd1; flipud(csf_usd1)];
fill(x1, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
plot(x, cdata_36hours1_global1, 'LineWidth', 2.0, 'Color', colormap_jet(1,:), 'DisplayName', sprintf('H1 (NRMSE: %.3f)', wrmse1_global1));
plot(x, cdata_36hours1_global2, 'LineWidth', 2.0, 'Color', colormap_jet(2,:), 'DisplayName', sprintf('H2 (NRMSE: %.3f)', wrmse1_global2));
plot(x, cdata_36hours1_global3, 'LineWidth', 2.0, 'Color', colormap_jet(3,:), 'DisplayName', sprintf('H3 (NRMSE: %.3f)', wrmse1_global3));
plot(x, cdata_36hours1_global4, 'LineWidth', 2.0, 'Color', colormap_jet(4,:), 'DisplayName', sprintf('H4 (NRMSE: %.3f)', wrmse1_global4));
plot(x, cdata_36hours1_global5, 'LineWidth', 2.0, 'Color', colormap_jet(5,:), 'DisplayName', sprintf('H5 (NRMSE: %.3f)', wrmse1_global5));
plot(time_exp1, csf_conc_exp1, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
model_points_global1 = interp1(x, cdata_36hours1_global1, time_exp1);
model_points_global2 = interp1(x, cdata_36hours1_global2, time_exp1);
model_points_global3 = interp1(x, cdata_36hours1_global3, time_exp1);
model_points_global4 = interp1(x, cdata_36hours1_global4, time_exp1);
model_points_global5 = interp1(x, cdata_36hours1_global5, time_exp1);
errorbar(time_exp1, csf_conc_exp1, csf_std1, 'k--', 'LineWidth', 2.0, 'HandleVisibility','off', 'CapSize', 6, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'Marker', 'o');
scatter(time_exp1, csf_conc_exp1, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
scatter(time_exp1, model_points_global1, 80, colormap_jet(1,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global2, 80, colormap_jet(2,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global3, 80, colormap_jet(3,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global4, 80, colormap_jet(4,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
scatter(time_exp1, model_points_global5, 80, colormap_jet(5,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
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
