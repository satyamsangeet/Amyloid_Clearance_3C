% Base params
p0 = struct('r_bc',0.038,'r_bp',0.014,'r_cp',0.00537,'r_p',0.427, 'A',16.203,'sigma_A',0.772,'sigma_bc',1.131, 'sigma_bp',1.768,'sigma_cp',6.100,'sigma_p',4.253);

% Model
function dydt = model(t, y, p)
    sw = (mod(t,24) >= 8 && mod(t,24) < 24);
    sl = 1 - sw;
    dydt = zeros(3,1);
    dydt(1) = (p.A*sw + p.sigma_A*p.A*sl) - (p.r_bc*sw + p.sigma_bc*p.r_bc*sl + p.r_bp*sw + p.sigma_bp*p.r_bp*sl)*y(1);
    dydt(2) = (p.r_bc*sw + p.sigma_bc*p.r_bc*sl)*y(1) - (p.r_cp*sw + p.sigma_cp*p.r_cp*sl)*y(2);
    dydt(3) = (p.r_bp*sw + p.sigma_bp*p.r_bp*sl)*y(1) + (p.r_cp*sw + p.sigma_cp*p.r_cp*sl)*y(2) - (p.r_p*sw + p.sigma_p*p.r_p*sl)*y(3);
end

function [t,w] = euler(F, T, y0, dt)
    t = (T(1):dt:T(2))';
    w = zeros(length(t),length(y0));
    w(1,:) = y0;
    for k = 1:length(t)-1
        w(k+1,:) = w(k,:) + F(t(k),w(k,:))'*dt;
    end
end

scenarios_r = {
    p0,                                        'Baseline';
    setfield(p0,'r_bc',  p0.r_bc*1.5),        'r_{bc} ×1.5 (Glymphatic +50%)';
    setfield(p0,'r_bc',  p0.r_bc*0.5),        'r_{bc} ×0.5 (Glymphatic -50%)';
    setfield(p0,'r_bp',  p0.r_bp*1.5),        'r_{bp} ×1.5 (BBB +50%)';
    setfield(p0,'r_bp',  p0.r_bp*0.5),        'r_{bp} ×0.5 (BBB -50%)';
    setfield(p0,'r_p',   p0.r_p*1.5),         'r_p ×1.5 (Peripheral +50%)';
    setfield(p0,'r_p',   p0.r_p*0.5),         'r_p ×0.5 (Peripheral -50%)';
};

scenarios_s = {
    p0,                                               'Baseline';
    setfield(p0,'sigma_bc', p0.sigma_bc*1.5),        '\sigma_{bc} ×1.5 (Glymphatic sleep +50%)';
    setfield(p0,'sigma_bc', p0.sigma_bc*0.5),        '\sigma_{bc} ×0.5 (Glymphatic sleep -50%)';
    setfield(p0,'sigma_bp', p0.sigma_bp*1.5),        '\sigma_{bp} ×1.5 (BBB sleep +50%)';
    setfield(p0,'sigma_bp', p0.sigma_bp*0.5),        '\sigma_{bp} ×0.5 (BBB sleep -50%)';
    setfield(p0,'sigma_p',  p0.sigma_p*1.5),         '\sigma_p ×1.5 (Peripheral sleep +50%)';
    setfield(p0,'sigma_p',  p0.sigma_p*0.5),         '\sigma_p ×0.5 (Peripheral sleep -50%)';
};

T = [0, 24*100]; y0 = [0,600,20]; dt = 0.01;
x_range  = [2336, 2384];
t_ticks  = 2336:4:2384;
t_labels = 0:4:48;
sleep_x  = [2352 2360 2376 2384];

titles = {'Brain ISF (y_1)','CSF (y_2)','Plasma (y_3)'};
comps  = [1, 2, 3];
t_sim = (T(1):dt:T(2))';

sols_r = cell(size(scenarios_r,1),1);
for s = 1:size(scenarios_r,1)
    p = scenarios_r{s,1};
    [~,w] = euler(@(t,y) model(t,y,p), T, y0, dt);
    sols_r{s} = w;
end

sols_s = cell(size(scenarios_s,1),1);
for s = 1:size(scenarios_s,1)
    p = scenarios_s{s,1};
    [~,w] = euler(@(t,y) model(t,y,p), T, y0, dt);
    sols_s{s} = w;
end

colors = lines(7);

function add_sleep_patches(sleep_x)
    yl = ylim;
    for i = 1:2:length(sleep_x)-1
        patch([sleep_x(i) sleep_x(i+1) sleep_x(i+1) sleep_x(i)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.7 0.9 1],'FaceAlpha',0.3,'EdgeColor','none', ...
              'HandleVisibility','off');
    end
end

function plot_3comp(t_sim, sols, idx_plot, scenarios, comps, titles, x_range, t_ticks, t_labels, sleep_x, colors, fig_title)
    figure('Position',[100 100 900 800]);
    for c = 1:3
        subplot(3,1,c); hold on;
        for k = 1:length(idx_plot)
            s = idx_plot(k);
            plot(t_sim, sols{s}(:,comps(c)), ...
                 'Color', colors(k,:), 'LineWidth', 2, ...
                 'DisplayName', scenarios{s,2});
        end
        xlim(x_range);
        add_sleep_patches(sleep_x);
        xticks(t_ticks); xticklabels(t_labels);
        xlabel('Time (h)','FontWeight','bold');
        ylabel('A\beta42 (pg/ml)','FontWeight','bold');
        title(titles{c},'FontWeight','bold');
        legend('Location','best'); grid on; hold off;
    end
    sgtitle(fig_title,'FontWeight','bold');
end

function plot_ratio(t_sim, sols, idx_plot, scenarios, x_range, t_ticks, t_labels, sleep_x, colors, fig_title)
    figure('Position',[100 100 900 400]);
    hold on;
    for k = 1:length(idx_plot)
        s = idx_plot(k);
        % Guard against near-zero plasma values
        ratio = sols{s}(:,1) ./ max(sols{s}(:,3), 1e-12);
        plot(t_sim, ratio, ...
             'Color', colors(k,:), 'LineWidth', 2, ...
             'DisplayName', scenarios{s,2});
    end
    xlim(x_range);
    add_sleep_patches(sleep_x);
    xticks(t_ticks); xticklabels(t_labels);
    xlabel('Time (h)','FontWeight','bold');
    ylabel('y_1 / y_3  (Brain ISF / Plasma)','FontWeight','bold');
    legend('Location','best'); grid on; hold off;
    sgtitle(fig_title,'FontWeight','bold');
end

% Plot 1: r_bc
plot_3comp(t_sim, sols_r, [1,2,3], scenarios_r, comps, titles, x_range, t_ticks, t_labels, sleep_x, colors, 'Effect of Glymphatic Clearance (r_{bc}) on All Compartments');

% Plot 2: r_p
plot_3comp(t_sim, sols_r, [1,6,7], scenarios_r, comps, titles, x_range, t_ticks, t_labels, sleep_x, colors, 'Effect of Systemic Clearance (r_p) on All Compartments');

% Plot 3: r_bp
plot_3comp(t_sim, sols_r, [1,4,5], scenarios_r, comps, titles, x_range, t_ticks, t_labels, sleep_x, colors, 'Effect of BBB Transport (r_{bp}) on All Compartments');

% Plot 4: sigma_bc
plot_3comp(t_sim, sols_s, [1,2,3], scenarios_s, comps, titles, x_range, t_ticks, t_labels, sleep_x, colors, 'Effect of Sleep-Phase Glymphatic Scaling (\sigma_{bc}) on All Compartments');

% Plot 5: sigma_bp
plot_3comp(t_sim, sols_s, [1,4,5], scenarios_s, comps, titles, x_range, t_ticks, t_labels, sleep_x, colors, 'Effect of Sleep-Phase BBB Scaling (\sigma_{bp}) on All Compartments');

% Plot 6: sigma_p
plot_3comp(t_sim, sols_s, [1,6,7], scenarios_s, comps, titles, x_range, t_ticks, t_labels, sleep_x, colors, 'Effect of Sleep-Phase Peripheral Scaling (\sigma_p) on All Compartments');

% Plot 7: ratio — r_bc
plot_ratio(t_sim, sols_r, [1,2,3], scenarios_r, x_range, t_ticks, t_labels, sleep_x, colors, 'Brain ISF / Plasma Ratio — Effect of r_{bc} (Glymphatic Clearance)');

% Plot 8: ratio — r_bp
plot_ratio(t_sim, sols_r, [1,4,5], scenarios_r, x_range, t_ticks, t_labels, sleep_x, colors, 'Brain ISF / Plasma Ratio — Effect of r_{bp} (BBB Transport)');

% Plot 9: ratio — r_p
plot_ratio(t_sim, sols_r, [1,6,7], scenarios_r, x_range, t_ticks, t_labels, sleep_x, colors, 'Brain ISF / Plasma Ratio — Effect of r_p (Systemic Clearance)');

% Plot 10: ratio — sigma_bc
plot_ratio(t_sim, sols_s, [1,2,3], scenarios_s, x_range, t_ticks, t_labels, sleep_x, colors, 'Brain ISF / Plasma Ratio — Effect of \sigma_{bc} (Glymphatic Sleep Scaling)');

% Plot 11: ratio — sigma_bp
plot_ratio(t_sim, sols_s, [1,4,5], scenarios_s, x_range, t_ticks, t_labels, sleep_x, colors, 'Brain ISF / Plasma Ratio — Effect of \sigma_{bp} (BBB Sleep Scaling)');

% Plot 12: ratio — sigma_p
plot_ratio(t_sim, sols_s, [1,6,7], scenarios_s, x_range, t_ticks, t_labels, sleep_x, colors, 'Brain ISF / Plasma Ratio — Effect of \sigma_p (Peripheral Sleep Scaling)');

figure('Position',[100 100 1100 500]);
[~,i16] = min(abs(t_sim - 2352));

labels_r = {'r_{bc}×1.5','r_{bc}×0.5','r_{bp}×1.5','r_{bp}×0.5','r_p×1.5','r_p×0.5'};
labels_s = {'\sigma_{bc}×1.5','\sigma_{bc}×0.5','\sigma_{bp}×1.5','\sigma_{bp}×0.5',...
            '\sigma_p×1.5','\sigma_p×0.5'};

idx_r = [2,3,4,5,6,7];
idx_s = [2,3,4,5,6,7];

pct_r = zeros(6,3);
pct_s = zeros(6,3);

for k = 1:6
    for c = 1:3
        ref = sols_r{1}(i16,c);
        pct_r(k,c) = (sols_r{idx_r(k)}(i16,c) - ref) / ref * 100;
        pct_s(k,c) = (sols_s{idx_s(k)}(i16,c) - ref) / ref * 100;
    end
end

all_labels = [labels_r, labels_s];
all_pct    = [pct_r; pct_s];   % 12×3
n = size(all_pct,1);
x = 1:n;
bw = 0.25;
bar_colors = [0.2 0.4 0.8; 0.2 0.7 0.5; 0.8 0.2 0.2];

hold on;
for c = 1:3
    bar(x + (c-2)*bw, all_pct(:,c), bw, ...
        'FaceColor', bar_colors(c,:), ...
        'DisplayName', titles{c});
end
yline(0,'k-','LineWidth',1);
xline(6.5, '--k', 'LineWidth', 1, 'Label','r_* | \sigma_*', 'LabelHorizontalAlignment','center', 'HandleVisibility','off');
xticks(x); xticklabels(all_labels); xtickangle(30);
ylabel('% Change from Baseline','FontWeight','bold');
title('Pathway & Sleep-Scaling Sensitivity — All Compartments at t=16h','FontWeight','bold');
legend('Location','best'); grid on; hold off;

fprintf('%-22s %12s %12s %12s\n','Intervention','Brain ISF','CSF','Plasma');
fprintf('%s\n', repmat('-',1,60));
for k = 1:6
    fprintf('%-22s %+11.1f%% %+11.1f%% %+11.1f%%\n', ...
            labels_r{k}, pct_r(k,1), pct_r(k,2), pct_r(k,3));
end

fprintf('%-22s %12s %12s %12s\n','Intervention','Brain ISF','CSF','Plasma');
fprintf('%s\n', repmat('-',1,60));
for k = 1:6
    fprintf('%-22s %+11.1f%% %+11.1f%% %+11.1f%%\n', ...
            strrep(labels_s{k},'\',''), pct_s(k,1), pct_s(k,2), pct_s(k,3));
end

fprintf('\n=== SUMMARY: y1/y3 RATIO AT t=16h ===\n');
fprintf('%-22s %12s %12s\n','Intervention','Ratio (y1/y3)','vs Baseline');
fprintf('%s\n', repmat('-',1,48));

base_ratio = sols_r{1}(i16,1) / max(sols_r{1}(i16,3), 1e-12);
fprintf('Baseline              %12.4f %+11s\n', base_ratio, '—');

all_scenario_labels = [labels_r, labels_s];
all_idx_r  = [idx_r, ones(1,6)];
for k = 1:6
    r = sols_r{idx_r(k)}(i16,1) / max(sols_r{idx_r(k)}(i16,3), 1e-12);
    fprintf('%-22s %12.4f %+10.1f%%\n', labels_r{k}, r, (r-base_ratio)/base_ratio*100);
end
for k = 1:6
    r = sols_s{idx_s(k)}(i16,1) / max(sols_s{idx_s(k)}(i16,3), 1e-12);
    fprintf('%-22s %12.4f %+10.1f%%\n', strrep(labels_s{k},'\',''), r, (r-base_ratio)/base_ratio*100);
end
