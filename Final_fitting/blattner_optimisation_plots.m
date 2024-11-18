% general graphics, this will apply to any figure you open
% (groot is the default figure object).
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 0.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 8, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 8, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'on', ...
'DefaultAxesTickLength', [0.02 0.025]);

% Updated model function to optimize only sigma_bp (a), sigma_cp (b), and rbc (a12_wake)
function dydt_n = model_global_wrmse(t, y)
    A = 13;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_cp = 0.0030;
    sigma_cp = 9.99;
    r_bp = (r_bc*(1-133*r_cp))/(133*r_cp);
    sigma_bp = 9.95;
    r_p = 0.28;
    sigma_p = 1.90;

    % Switch between sleep and wake states
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_cp * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model_global_rmse(t, y)
    A = 13;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_cp = 0.0029;
    sigma_cp = 9.42;
    r_bp = (r_bc*(1-133*r_cp))/(133*r_cp);
    sigma_bp = 9.78;
    r_p = 0.28;
    sigma_p = 1.88;

    % Switch between sleep and wake states
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_cp * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model_ind_wrmse(t, y)
    A = 13;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_cp = 0.0042;
    sigma_cp = 4.99;
    r_bp = (r_bc*(1-133*r_cp))/(133*r_cp);
    sigma_bp = 1.01;
    r_p = 0.28;
    sigma_p = 2.99;

    % Switch between sleep and wake states
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_cp * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

function dydt_n = model_ind_rmse(t, y)
    A = 13;
    sigma_A = 0.8;
    r_bc = 1.5;
    sigma_bc = 2.5;
    r_cp = 0.002;
    sigma_cp = 9.98;
    r_bp = (r_bc*(1-133*r_cp))/(133*r_cp);
    sigma_bp = 1.05;
    r_p = 0.28;
    sigma_p = 5.44;

    % Switch between sleep and wake states
    sw_cycle = (mod(t, 24) >= 8 && mod(t, 24) < 24);
    
    % Define the ODE system
    dydt_n = zeros(3, 1);
    dydt_n(1) = A * sw_cycle + sigma_A * A * (1 - sw_cycle) - (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle) + r_bc * sw_cycle + sigma_cp * r_bc * (1 - sw_cycle)) * y(1);
    dydt_n(2) = (r_bc * sw_cycle + sigma_bc * r_bc * (1 - sw_cycle)) * y(1) - (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2);
    dydt_n(3) = (r_bp * sw_cycle + sigma_bp * r_bp * (1 - sw_cycle)) * y(1) + (r_cp * sw_cycle + sigma_cp * r_cp * (1 - sw_cycle)) * y(2) - (r_p * sw_cycle + sigma_p * r_p * (1 - sw_cycle)) * y(3);
end

% Defining Euler-Maruyama Method
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

[t_100days_gwrmse, sol_100days_gwrmse] = euler(@(t,y) model_global_wrmse(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_grmse, sol_100days_grmse] = euler(@(t,y) model_global_rmse(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_iwrmse, sol_100days_iwrmse] = euler(@(t,y) model_ind_wrmse(t,y), [0, 24*100], [0,600,15.5], 0.01);
[t_100days_irmse, sol_100days_irmse] = euler(@(t,y) model_ind_rmse(t,y), [0, 24*100], [0,600,15.5], 0.01);

cdata_12hours1 = sol_100days_gwrmse(233600:234800, 2);
cdata_12hours2 = sol_100days_grmse(233600:234800, 2);
cdata_12hours3 = sol_100days_iwrmse(233600:234800, 2);
cdata_12hours4 = sol_100days_irmse(233600:234800, 2);
cdata_36hours1 = sol_100days_gwrmse(233600:237200, 2);
cdata_36hours2 = sol_100days_grmse(233600:237200, 2);
cdata_36hours3 = sol_100days_iwrmse(233600:237200, 2);
cdata_36hours4 = sol_100days_irmse(233600:237200, 2);
csf_mean1 = mean(cdata_12hours1);
csf_mean2 = mean(cdata_12hours2);
csf_mean3 = mean(cdata_12hours3);
csf_mean4 = mean(cdata_12hours4);
normalise_C1 = (cdata_36hours1 / csf_mean1)*100;
normalise_C2 = (cdata_36hours2 / csf_mean2)*100;
normalise_C3 = (cdata_36hours3 / csf_mean3)*100;
normalise_C4 = (cdata_36hours4 / csf_mean4)*100;

experimental_csf = readtable('/home/satyam/Documents/combined_fitting_data/final_fitting/blattner2020_csf42.csv');
experimental_csf.Time = experimental_csf.Time;
time_exp = experimental_csf.Time;
csf_conc_exp = experimental_csf.Conc;
lsd_csf = experimental_csf.LSD;
usd_csf = experimental_csf.USD;

x = 1:3601;
ticks = 101:100:3601;
x = 1:3601;
disp(length(x));

% Define new ticks and labels
new_ticks = 1:100:3600+100;
new_labels = 1:length(new_ticks);

% Ensure ticks + 1 does not exceed the length of normalise_C
valid_indices = ticks + 1 <= length(normalise_C1);
valid_ticks = ticks(valid_indices);

% Plot comparison CSF
figure();
x2 = [time_exp; flipud(time_exp)];
inBetween1 = [lsd_csf; flipud(usd_csf)];
fill(x2, inBetween1, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'k', 'HandleVisibility', 'off');
hold on;
%plot(x, normalise_C1, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Global WRMSE');
%plot(x, normalise_C2, 'g--', 'LineWidth', 2.0, 'DisplayName', 'Global RMSE');
plot(x, normalise_C3, 'b--', 'LineWidth', 2.0, 'DisplayName', 'Ind WRMSE');
%plot(x, normalise_C4, 'm--', 'LineWidth', 2.0, 'DisplayName', 'Ind RMSE');
plot(time_exp, csf_conc_exp, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Blattner2020 - CSF');
scatter(time_exp, csf_conc_exp, 'ko', 'MarkerFaceColor','k', 'HandleVisibility', 'off');
xline(1600, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xline(2400, 'k--', 'LineWidth', 1.0, 'Alpha', 0.3, 'HandleVisibility', 'off');
xlim([2336, 2386]);
xticks(0:200:3600);
xticklabels(0:2:50);
legend('show');
xlabel('Time (hr)', 'FontWeight', 'bold');
ylabel('Percetage Baseline CSF', 'FontWeight', 'bold');
xlim([0, 3500]);
hold off;