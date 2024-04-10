
%% Two Compartment Model for Brain clearance
%% Compartment 1 = Blood (B)
%% Compartment 2 = CSF (C)

%% Model Parameters
a12 = 0.3;                              % Amyloid transfer from Blood to CSF
a21_wake = 0.1;                         % Amyloid transfer from CSF to blood during wake state
a21_sleep = 0.5;             % Amyloid transfer from CSF to blood during sleep state
A_wake = 69;                            % Amyloid production during wake
A_sleep = 1;                            % Amyloid production during sleep
k = 2;                                  % Amyloid clearance from blood

% Time span
t_start = 0.0;
t_end = 24*5;                           % Time in hours (24 hours X 5 days)
dt = 0.001;

initial_conditions = [0;0];             % Initial conditions for Blood and CSF compartment respectively

% Defining the model
ode_system1 = @(t,y1) model1(t, y1, a12, a21_wake, a21_sleep, A_wake, A_sleep, k);
ode_system2 = @(t,y2) model2(t, y2, a21_wake, a21_sleep, A_wake, A_sleep, k);
options = odeset('MaxStep', dt);
[t,y1] = ode45(ode_system1, t_start:dt:t_end, initial_conditions, options);
[t,y2] = ode45(ode_system2, t_start:dt:t_end, initial_conditions, options);

%Plasma to CSF ratio
ratio1 = y1(:,1) ./ y1(:,2);
ratio2 = y2(:,1) ./ y2(:,2);

% Plotting
figure;
subplot(3, 1, 1);
hold on;
plot(t, y1(:, 1), 'LineWidth', 3, 'DisplayName', 'Plasma (with a12)');
plot(t, y2(:, 1), 'LineWidth', 3, 'DisplayName', 'Plasma (without a12)');
ylabel('Concentration (pg/ml)');
set(gca, 'XTick', [0 24 48 72 96 120]);
set(gca, 'XTickLabel', {'0', '24', '48', '72', '96', '120'});
legend('show');

subplot(3, 1, 2);
hold on;
plot(t, y1(:, 2), 'LineWidth', 3, 'DisplayName', 'CSF (With a12)');
plot(t, y2(:, 2), 'LineWidth', 3, 'DisplayName', 'CSF (without a12)');
xlabel('Time (hrs)');
ylabel('Concentration (pg/ml)');
set(gca, 'XTick', [0 24 48 72 96 120]);
set(gca, 'XTickLabel', {'0', '24', '48', '72', '96', '120'});
legend('show');

subplot(3, 1, 3);
hold on;
plot(t, ratio1, 'LineWidth', 3, 'DisplayName', 'Plasma/CSF Ratio with a12');
plot(t, ratio2, 'LineWidth', 3, 'DisplayName', 'Plasma/CSF Ratio without a12');
xlabel('Time (hrs)');
ylabel('Ratio');
set(gca, 'XTick', [0 24 48 72 96 120]);
set(gca, 'XTickLabel', {'0', '24', '48', '72', '96', '120'});
legend('show');

% Steady state solution with a12
%for wake state
Bw12 = A_wake/k;
Cw12 = ((a12 + k) * A_wake)/(a21_wake * k);

%for sleep state
Bs12 = A_sleep/k;
Cs12 = ((a12 + k)* A_sleep)/(a21_sleep * k);

%Steady state solution without a12
Bw = A_wake/k;
Cw = A_wake/a21_wake;

Bs = A_sleep/k;
Cs = A_sleep/a21_sleep;

fprintf('The steady state value with parameter a12:: \n');
fprintf('For Wake state:: \n');
fprintf('Blood Compartment : %f', Bw12);
fprintf('\nCSF Compartment : %f', Cw12);
fprintf('\nFor Sleep state:: \n');
fprintf('Blood Compartment : %f', Bs12);
fprintf('\nCSF Compartment : %f', Cs12);

fprintf('\nThe steady state value without parameter a12:: \n');
fprintf('For Wake state:: \n');
fprintf('Blood Compartment : %f', Bw);
fprintf('\nCSF Compartment : %f', Cw);
fprintf('\nFor Sleep state:: \n');
fprintf('Blood Compartment : %f', Bs);
fprintf('\nCSF Compartment : %f', Cs);

% Defining the ODEs
function dydt1 = model1(t, y1, a12, a21_wake, a21_sleep, A_wake, A_sleep, k)
    sw_cycle = (mod(t,24) >= 8 && mod(t,24) < 24);
    dydt1 = zeros(2, 1);
    dydt1(1) = (a21_wake * sw_cycle + a21_sleep * (1 - sw_cycle)) * y1(2) - (a12 + k) * y1(1);
    dydt1(2) = a12 * y1(1) + A_wake * sw_cycle + A_sleep * (1 - sw_cycle) - (a21_wake * sw_cycle + a21_sleep * (1 - sw_cycle)) * y1(2);
end

function dydt2 = model2(t, y2, a21_wake, a21_sleep, A_wake, A_sleep, k)
    sw_cycle = (mod(t,24) >= 8 && mod(t,24) < 24);
    dydt2 = zeros(2, 1);
    dydt2(1) = (a21_wake * sw_cycle + a21_sleep * (1 - sw_cycle)) * y2(2) - k * y2(1);
    dydt2(2) = A_wake * sw_cycle + A_sleep * (1 - sw_cycle) - (a21_wake * sw_cycle + a21_sleep * (1 - sw_cycle)) * y2(2);
end
