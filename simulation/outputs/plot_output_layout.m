clc
close all
clear

%% Configuration
sensor_dist = 0.5;
sensor_str = num2str(sensor_dist);
%sensor_str = 'Delta1';

% Figure size
fig_x = 50;
fig_y = 50;
fig_w = 560;
fig_h = 400;

% Y-axis limits
axle_MAE_ylim  = [0 20];
axle_RMSE_ylim = [0 25];
axle_MRE_ylim  = [0 55];
gvw_MAE_ylim   = [0 15];
gvw_RMSE_ylim  = [0 18];
gvw_MRE_ylim   = [0 35];

% Number of sensors
sensor_counts = [2, 4, 6, 8, 10, 12, 14, 16];
markers = {'-o', '-s', '-^', '-d', '-p', '-h', '-x', '-+'};

% Load data from CSV files
axle_data = cell(length(sensor_counts), 1);
gvw_data = cell(length(sensor_counts), 1);
for i = 1:length(sensor_counts)
    sensors = sensor_counts(i);
    axle_data{i} = readtable(sprintf('axle_output_s%d_d%s.csv', sensors, sensor_str));
    gvw_data{i} = readtable(sprintf('gvw_output_s%d_d%s.csv', sensors, sensor_str));
end

%% Axle maximum relative error
figure;
hold on;
%title(sprintf('Axle maximum relative error, d = %sm',sensor_str));
xlabel('Number of Sensors');
ylabel('Maximum Relative Error (%)');
grid on;
grid minor;
ylim(axle_MRE_ylim);
set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    max_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = axle_data{i};
        max_mv_values(i) = df.max_mv(df.speed == speed);
    end
    p = plot(sensor_counts, max_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
saveas(gcf, sprintf('figures_layout/axle_max_er_d%s.png',sensor_str));

%% GVW maximum relative error
figure;
hold on;
%title(sprintf('GVW maximum relative error, d = %sm',sensor_str));
xlabel('Number of Sensors');
ylabel('Maximum Relative Error (%)');
grid on;
grid minor;
ylim(gvw_MRE_ylim);
set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    max_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = gvw_data{i}; 
        max_mv_values(i) = df.max_mv(df.speed == speed);
    end
    p = plot(sensor_counts, max_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
saveas(gcf, sprintf('figures_layout/gvw_max_er_d%s.png',sensor_str));

%% Axle mean absolute error
figure;
hold on;
%title(sprintf('Axle MAE, d = %sm',sensor_str));
xlabel('Number of Sensors');
ylabel('MAE (%)');
grid on;
grid minor;
ylim(axle_MAE_ylim);
set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    mae_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = axle_data{i};
        mae_mv_values(i) = df.mae_mv(df.speed == speed);
    end
    p = plot(sensor_counts, mae_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
saveas(gcf, sprintf('figures_layout/axle_mae_d%s.png',sensor_str));

%% GVW mean absolute error
figure;
hold on;
%title(sprintf('GVW MAE, d = %sm',sensor_str));
xlabel('Number of Sensors');
ylabel('MAE (%)');
grid on;
grid minor;
ylim(gvw_MAE_ylim);
set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    mae_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = gvw_data{i};
        mae_mv_values(i) = df.mae_mv(df.speed == speed);
    end
    p = plot(sensor_counts, mae_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
saveas(gcf, sprintf('figures_layout/gvw_mae_d%s.png',sensor_str));

%% Axle root mean squared error
figure;
hold on;
%title(sprintf('Axle RMSE, d = %sm',sensor_str));
xlabel('Number of Sensors');
ylabel('RMSE (%)');
grid on;
grid minor;
ylim(axle_RMSE_ylim);
set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    rmse_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = axle_data{i};
        rmse_mv_values(i) = df.rmse_mv(df.speed == speed);
    end
    p = plot(sensor_counts, rmse_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
saveas(gcf, sprintf('figures_layout/axle_rmse_d%s.png',sensor_str));

%% GVW root mean squared error
figure;
hold on;
%title(sprintf('GVW RMSE, d = %sm',sensor_str));
xlabel('Number of Sensors');
ylabel('RMSE (%)');
grid on;
grid minor;
ylim(gvw_RMSE_ylim);
set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    rmse_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = gvw_data{i};
        rmse_mv_values(i) = df.rmse_mv(df.speed == speed);
    end
    p = plot(sensor_counts, rmse_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
saveas(gcf, sprintf('figures_layout/gvw_rmse_d%s.png',sensor_str));

return;

%% Axle standard deviation
% figure;
% hold on;
% title(sprintf('Axle standard deviation, d = %sm',sensor_str));
% xlabel('Number of Sensors');
% ylabel('Standard deviation');
% grid on;
% 
% unique_speeds = unique(axle_data{1}.speed);
% for speed_idx = 1:length(unique_speeds)
%     speed = unique_speeds(speed_idx);
%     std_mv_values = zeros(length(sensor_counts), 1);
%     for i = 1:length(sensor_counts)
%         df = axle_data{i};
%         std_mv_values(i) = df.std_mv(df.speed == speed);
%     end
%     p = plot(sensor_counts, std_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
%     p.MarkerFaceColor = p.Color;
% end
% legend;
% hold off;
% saveas(gcf, sprintf('figures_layout/axle_std_d%s.png',sensor_str));

%% GVW standard deviation
% figure;
% hold on;
% title(sprintf('GVW standard deviation, d = %sm',sensor_str));
% xlabel('Number of Sensors');
% ylabel('Standard deviation');
% grid on;
% 
% unique_speeds = unique(axle_data{1}.speed);
% for speed_idx = 1:length(unique_speeds)
%     speed = unique_speeds(speed_idx);
%     std_mv_values = zeros(length(sensor_counts), 1);
%     for i = 1:length(sensor_counts)
%         df = gvw_data{i};
%         std_mv_values(i) = df.std_mv(df.speed == speed);
%     end
%     p = plot(sensor_counts, std_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
%     p.MarkerFaceColor = p.Color;
% end
% legend;
% hold off;
% saveas(gcf, sprintf('figures_layout/gvw_std_d%s.png',sensor_str));
