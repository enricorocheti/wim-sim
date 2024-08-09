clc;
close all;
clear;

%% Configuration
sensor_dist = '3';
sensor_counts = [2, 4, 6, 8, 10, 12, 14, 16];
markers = {'-o', '-s', '-^', '-d', '-p', '-h', '-x', '-+'};
estimators = {'mv', 'MLE', 'pchip', 'makima', 'spline'};
speed_labels = {'20', '40', '60', '80', '100'};

% Figure size
fig_x = 50;
fig_y = 50;
fig_w = 560;
fig_h = 400;
fontSize = 17;
fontSize2 = 12;

% Y-axis limits
axle_MRE_ylim  = [0 50];
axle_MAE_ylim  = [0 16];
axle_RMSE_ylim = [0 20];
gvw_MRE_ylim   = [0 28];
gvw_MAE_ylim   = [0 12];
gvw_RMSE_ylim  = [0 14];

% Load data from CSV files
axle_data = cell(length(sensor_counts), 1);
gvw_data = cell(length(sensor_counts), 1);
for i = 1:length(sensor_counts)
    sensors = sensor_counts(i);
    axle_data{i} = readtable(sprintf('axle_output_s%d_d%s.csv', sensors, sensor_dist));
    gvw_data{i} = readtable(sprintf('gvw_output_s%d_d%s.csv', sensors, sensor_dist));
end

unique_speeds = unique(axle_data{1}.speed);

for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    speed_label = speed_labels{speed_idx};

    %% Axle maximum relative error
    figure;
    hold on;
    %title(sprintf('Axle maximum relative error, Speed = %d km/h, d = %s', speed, sensor_dist));
    set(gca, 'FontSize', fontSize2);
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('Maximum Relative Error (%)', 'FontSize', fontSize);
    grid on;
    grid minor;
    %ylim(axle_MRE_ylim);
    set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        max_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            max_values(i) = df.(sprintf('max_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, max_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/axle_max_er_v%s_%s.png', speed_label, sensor_dist));

    %% GVW maximum relative error
    figure;
    hold on;
    set(gca, 'FontSize', fontSize2);
    %title(sprintf('GVW maximum relative error, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('Maximum Relative Error (%)', 'FontSize', fontSize);
    grid on;
    grid minor;
    %ylim(gvw_MRE_ylim);
    set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        max_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            max_values(i) = df.(sprintf('max_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, max_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/gvw_max_er_v%s_%s.png', speed_label, sensor_dist));

%     %% Axle mean absolute error
%     figure;
%     hold on;
%     set(gca, 'FontSize', fontSize2);
%     %title(sprintf('Axle MAE, Speed = %d km/h, d = %s', speed, sensor_dist));
%     xlabel('Number of Sensors');
%     ylabel('MAE (%)');
%     grid on;
%     grid minor;
%     %ylim(axle_MAE_ylim);
%     set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);
% 
%     for estimator_idx = 1:length(estimators)
%         estimator = estimators{estimator_idx};
%         mae_values = zeros(length(sensor_counts), 1);
%         for i = 1:length(sensor_counts)
%             df = axle_data{i};
%             mae_values(i) = df.(sprintf('mae_%s', estimator))(df.speed == speed);
%         end
%         p = plot(sensor_counts, mae_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
%         p.MarkerFaceColor = p.Color;
%     end
%     legend;
%     hold off;
%     saveas(gcf, sprintf('figures_estimators/group_speed/axle_mae_v%s_%s.png', speed_label, sensor_dist));
% 
%     %% GVW mean absolute error
%     figure;
%     hold on;
%     set(gca, 'FontSize', fontSize2);
%     %title(sprintf('GVW MAE, Speed = %d km/h, d = %s', speed, sensor_dist));
%     xlabel('Number of Sensors');
%     ylabel('MAE (%)');
%     grid on;
%     grid minor;
%     %ylim(gvw_MAE_ylim);
%     set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);
% 
%     for estimator_idx = 1:length(estimators)
%         estimator = estimators{estimator_idx};
%         mae_values = zeros(length(sensor_counts), 1);
%         for i = 1:length(sensor_counts)
%             df = gvw_data{i};
%             mae_values(i) = df.(sprintf('mae_%s', estimator))(df.speed == speed);
%         end
%         p = plot(sensor_counts, mae_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
%         p.MarkerFaceColor = p.Color;
%     end
%     legend;
%     hold off;
%     saveas(gcf, sprintf('figures_estimators/group_speed/gvw_mae_v%s_%s.png', speed_label, sensor_dist));
    
    %% Axle root mean squared error
    figure;
    hold on;
    set(gca, 'FontSize', fontSize2);
    %title(sprintf('Axle RMSE, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('RMSE (%)', 'FontSize', fontSize);
    grid on;
    grid minor;
    %ylim(axle_RMSE_ylim);
    set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        rmse_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            rmse_values(i) = df.(sprintf('rmse_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, rmse_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/axle_rmse_v%s_%s.png', speed_label, sensor_dist));

    %% GVW root mean squared error
    figure;
    hold on;
    set(gca, 'FontSize', fontSize2);
    %title(sprintf('GVW RMSE, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('RMSE (%)', 'FontSize', fontSize);
    grid on;
    grid minor;
    %ylim(gvw_RMSE_ylim);
    set(gcf,'position',[fig_x,fig_y,fig_w,fig_h]);

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        rmse_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            rmse_values(i) = df.(sprintf('rmse_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, rmse_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/gvw_rmse_v%s_%s.png', speed_label, sensor_dist));

    %% Axle standard deviation
%     figure;
%     hold on;
%     %title(sprintf('Axle standard deviation, Speed = %d km/h, d = %s', speed, sensor_dist));
%     xlabel('Number of Sensors');
%     ylabel('Standard deviation');
%     grid on;
% 
%     for estimator_idx = 1:length(estimators)
%         estimator = estimators{estimator_idx};
%         std_values = zeros(length(sensor_counts), 1);
%         for i = 1:length(sensor_counts)
%             df = axle_data{i};
%             std_values(i) = df.(sprintf('std_%s', estimator))(df.speed == speed);
%         end
%         p = plot(sensor_counts, std_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
%         p.MarkerFaceColor = p.Color;
%     end
%     legend;
%     hold off;
%     saveas(gcf, sprintf('figures_estimators/group_speed/axle_std_v%s_%s.png', speed_label, sensor_dist));

    %% GVW standard deviation
%     figure;
%     hold on;
%     %title(sprintf('GVW standard deviation, Speed = %d km/h, d = %s', speed, sensor_dist));
%     xlabel('Number of Sensors');
%     ylabel('Standard deviation');
%     grid on;
% 
%     for estimator_idx = 1:length(estimators)
%         estimator = estimators{estimator_idx};
%         std_values = zeros(length(sensor_counts), 1);
%         for i = 1:length(sensor_counts)
%             df = gvw_data{i};
%             std_values(i) = df.(sprintf('std_%s', estimator))(df.speed == speed);
%         end
%         p = plot(sensor_counts, std_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
%         p.MarkerFaceColor = p.Color;
%     end
%     legend;
%     hold off;
%     saveas(gcf, sprintf('figures_estimators/group_speed/gvw_std_v%s_%s.png', speed_label, sensor_dist));
end