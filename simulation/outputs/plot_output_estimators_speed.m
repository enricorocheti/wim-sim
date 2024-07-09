clc;
close all;
clear;

%% Configuration
sensor_dist = 'Delta1';
sensor_counts = [2, 4, 6, 8, 10, 12, 14, 16];
markers = {'-o', '-s', '-^', '-d', '-p', '-h', '-x', '-+'};
estimators = {'mv', 'MLE', 'pchip', 'makima', 'spline'};
speed_labels = {'20', '40', '60', '80', '100'};

axle_data = cell(length(sensor_counts), 1);
gvw_data = cell(length(sensor_counts), 1);

% Load data from CSV files
for i = 1:length(sensor_counts)
    sensors = sensor_counts(i);
    axle_data{i} = readtable(sprintf('axle_output_s%d_%s.csv', sensors, sensor_dist));
    gvw_data{i} = readtable(sprintf('gvw_output_s%d_%s.csv', sensors, sensor_dist));
end

unique_speeds = unique(axle_data{1}.speed);

for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    speed_label = speed_labels{speed_idx};

    %% Axle maximum relative error
    figure;
    hold on;
    title(sprintf('Axle maximum relative error, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Er (%)');
    grid on;

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
    title(sprintf('GVW maximum relative error, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Er (%)');
    grid on;

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

    %% Axle absolute mean error
    figure;
    hold on;
    title(sprintf('Axle absolute mean error, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Mean error (%)');
    grid on;

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        mean_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            mean_values(i) = df.(sprintf('mean_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, mean_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/axle_abs_mean_v%s_%s.png', speed_label, sensor_dist));

    %% GVW absolute mean error
    figure;
    hold on;
    title(sprintf('GVW absolute mean error, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Mean error (%)');
    grid on;

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        mean_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            mean_values(i) = df.(sprintf('mean_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, mean_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/gvw_abs_mean_v%s_%s.png', speed_label, sensor_dist));

    %% Axle standard deviation
    figure;
    hold on;
    title(sprintf('Axle standard deviation, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Standard deviation');
    grid on;

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        std_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            std_values(i) = df.(sprintf('std_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, std_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/axle_std_v%s_%s.png', speed_label, sensor_dist));

    %% GVW standard deviation
    figure;
    hold on;
    title(sprintf('GVW standard deviation, Speed = %d km/h, d = %s', speed, sensor_dist));
    xlabel('Number of Sensors');
    ylabel('Standard deviation');
    grid on;

    for estimator_idx = 1:length(estimators)
        estimator = estimators{estimator_idx};
        std_values = zeros(length(sensor_counts), 1);
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            std_values(i) = df.(sprintf('std_%s', estimator))(df.speed == speed);
        end
        p = plot(sensor_counts, std_values, markers{estimator_idx}, 'DisplayName', sprintf('%s', estimator));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_estimators/group_speed/gvw_std_v%s_%s.png', speed_label, sensor_dist));
end