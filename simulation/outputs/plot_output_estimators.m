clc;
close all;
clear;

%% Configuration
sensor_dist = '1';
sensor_counts = [2, 4, 6, 8, 10, 12, 14, 16];
markers = {'-o', '-s', '-^', '-d', '-p', '-h', '-x', '-+'};
estimators = {'mv', 'MLE', 'pchip', 'makima', 'spline'};

% Figure size
fig_x = 50;
fig_y = 50;
fig_w = 700;
fig_h = 500;

% Y-axis limits
axle_MAE_ylim  = [0 20];
axle_RMSE_ylim = [0 25];
axle_MRE_ylim  = [0 55];
gvw_MAE_ylim   = [0 15];
gvw_RMSE_ylim  = [0 18];
gvw_MRE_ylim   = [0 35];

% Load data from CSV files
axle_data = cell(length(sensor_counts), 1);
gvw_data = cell(length(sensor_counts), 1);
for i = 1:length(sensor_counts)
    sensors = sensor_counts(i);
    axle_data{i} = readtable(sprintf('axle_output_s%d_d%s.csv', sensors, sensor_dist));
    gvw_data{i} = readtable(sprintf('gvw_output_s%d_d%s.csv', sensors, sensor_dist));
end

for estimator_idx = 1:length(estimators)
    estimator = estimators{estimator_idx};

    %% Axle maximum relative error
    figure;
    hold on;
    title(sprintf('Axle maximum relative error (%s), d = %s', estimator, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Er (%)');
    grid on;

    unique_speeds = unique(axle_data{1}.speed);
    for speed_idx = 1:length(unique_speeds)
        speed = unique_speeds(speed_idx);
        max_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            max_values(i) = df.(sprintf('max_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, max_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_estimator/axle_max_er_%s_%s.png', estimator, sensor_dist));

    %% GVW maximum relative error
    figure;
    hold on;
    title(sprintf('GVW maximum relative error (%s), d = %s', estimator, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Er (%)');
    grid on;

    unique_speeds = unique(gvw_data{1}.speed);
    for speed_idx = 1:length(unique_speeds)
        speed = unique_speeds(speed_idx);
        max_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            max_values(i) = df.(sprintf('max_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, max_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_estimator/gvw_max_er_%s_%s.png', estimator, sensor_dist));

    %% Axle absolute mean error
    figure;
    hold on;
    title(sprintf('Axle absolute mean error (%s), d = %s', estimator, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Mean error (%)');
    grid on;

    unique_speeds = unique(axle_data{1}.speed);
    for speed_idx = 1:length(unique_speeds)
        speed = unique_speeds(speed_idx);
        mean_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            mean_values(i) = df.(sprintf('mean_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, mean_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_estimator/axle_abs_mean_%s_%s.png', estimator, sensor_dist));

    %% GVW absolute mean error
    figure;
    hold on;
    title(sprintf('GVW absolute mean error (%s), d = %s', estimator, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Mean error (%)');
    grid on;

    unique_speeds = unique(gvw_data{1}.speed);
    for speed_idx = 1:length(unique_speeds)
        speed = unique_speeds(speed_idx);
        mean_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            mean_values(i) = df.(sprintf('mean_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, mean_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_estimator/gvw_abs_mean_%s_%s.png', estimator, sensor_dist));

    %% Axle standard deviation
    figure;
    hold on;
    title(sprintf('Axle standard deviation (%s), d = %s', estimator, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Standard deviation');
    grid on;

    unique_speeds = unique(axle_data{1}.speed);
    for speed_idx = 1:length(unique_speeds)
        speed = unique_speeds(speed_idx);
        std_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            std_values(i) = df.(sprintf('std_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, std_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_estimator/axle_std_%s_%s.png', estimator, sensor_dist));

    %% GVW standard deviation
    figure;
    hold on;
    title(sprintf('GVW standard deviation (%s), d = %s', estimator, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Standard deviation');
    grid on;

    unique_speeds = unique(gvw_data{1}.speed);
    for speed_idx = 1:length(unique_speeds)
        speed = unique_speeds(speed_idx);
        std_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            std_values(i) = df.(sprintf('std_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, std_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_estimator/gvw_std_%s_%s.png', estimator, sensor_dist));
end