clc;
close all;
clear;

%% Configuration
sensor_dists = {'0_5', '1', '2', '3', '4', '5', 'Delta1', 'Delta2'};
sensor_counts = [2, 4, 6, 8, 10, 12, 14, 16];
markers = {'-o', '-s', '-^', '-d', '-p', '-h', '-x', '-+'};
speed_labels = {'20', '40', '60', '80', '100'};

% Figure size
fig_x = 50;
fig_y = 50;
fig_w = 540;
fig_h = 400;
fontSize = 17;
fontSize2 = 12;

data = struct();

% Load data from CSV files
for dist_idx = 1:length(sensor_dists)
    sensor_dist = strrep(sensor_dists{dist_idx}, '_', '.');
    axle_data = cell(length(sensor_counts), 1);
    gvw_data = cell(length(sensor_counts), 1);
    
    for i = 1:length(sensor_counts)
        sensors = sensor_counts(i);
        axle_data{i} = readtable(sprintf('axle_output_s%d_d%s.csv', sensors, sensor_dist));
        gvw_data{i} = readtable(sprintf('gvw_output_s%d_d%s.csv', sensors, sensor_dist));
    end
    
    data.(['d' sensor_dists{dist_idx}]).axle_data = axle_data;
    data.(['d' sensor_dists{dist_idx}]).gvw_data = gvw_data;
end

unique_speeds = unique(data.d0_5.axle_data{1}.speed);

for selected_speed = [20, 40, 60, 80, 100]
    if ~ismember(selected_speed, unique_speeds)
        error('Selected speed %d km/h is not available in the dataset.', selected_speed);
    end

    %% Axle maximum relative error
    figure('Position', [fig_x, fig_y, fig_w, fig_h]);
    hold on;
    set(gca, 'FontSize', fontSize2);
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('Maximum Relative Error (%)', 'FontSize', fontSize);
    grid on;
    grid minor;

    for dist_idx = 1:length(sensor_dists)
        sensor_dist = sensor_dists{dist_idx};
        axle_data = data.(['d' sensor_dist]).axle_data;
        max_values = zeros(length(sensor_counts), 1);
        
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            max_values(i) = df.max_mv(df.speed == selected_speed);
        end
        p = plot(sensor_counts, max_values, markers{dist_idx}, 'DisplayName', sprintf('d = %s', strrep(sensor_dist, '_', '.')));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_layout/axle_max_er_v%s.png', num2str(selected_speed)));

    %% GVW maximum relative error
    figure('Position', [fig_x, fig_y, fig_w, fig_h]);
    hold on;
    set(gca, 'FontSize', fontSize2);
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('Maximum Relative Error (%)', 'FontSize', fontSize);
    grid on;
    grid minor;

    for dist_idx = 1:length(sensor_dists)
        sensor_dist = sensor_dists{dist_idx};
        gvw_data = data.(['d' sensor_dist]).gvw_data;
        max_values = zeros(length(sensor_counts), 1);
        
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            max_values(i) = df.max_mv(df.speed == selected_speed);
        end
        p = plot(sensor_counts, max_values, markers{dist_idx}, 'DisplayName', sprintf('d = %s', strrep(sensor_dist, '_', '.')));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_layout/gvw_max_er_v%s.png', num2str(selected_speed)));

    %% Axle RMSE
    figure('Position', [fig_x, fig_y, fig_w, fig_h]);
    hold on;
    set(gca, 'FontSize', fontSize2);
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('RMSE (%)', 'FontSize', fontSize);
    grid on;
    grid minor;

    for dist_idx = 1:length(sensor_dists)
        sensor_dist = sensor_dists{dist_idx};
        axle_data = data.(['d' sensor_dist]).axle_data;
        rmse_values = zeros(length(sensor_counts), 1);
        
        for i = 1:length(sensor_counts)
            df = axle_data{i};
            rmse_values(i) = df.rmse_mv(df.speed == selected_speed);
        end
        p = plot(sensor_counts, rmse_values, markers{dist_idx}, 'DisplayName', sprintf('d = %s', strrep(sensor_dist, '_', '.')));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_layout/axle_rmse_v%s.png', num2str(selected_speed)));

    %% GVW RMSE
    figure('Position', [fig_x, fig_y, fig_w, fig_h]);
    hold on;
    set(gca, 'FontSize', fontSize2);
    xlabel('Number of Sensors', 'FontSize', fontSize);
    ylabel('RMSE (%)', 'FontSize', fontSize);
    grid on;
    grid minor;

    for dist_idx = 1:length(sensor_dists)
        sensor_dist = sensor_dists{dist_idx};
        gvw_data = data.(['d' sensor_dist]).gvw_data;
        rmse_values = zeros(length(sensor_counts), 1);
        
        for i = 1:length(sensor_counts)
            df = gvw_data{i};
            rmse_values(i) = df.rmse_mv(df.speed == selected_speed);
        end
        p = plot(sensor_counts, rmse_values, markers{dist_idx}, 'DisplayName', sprintf('d = %s', strrep(sensor_dist, '_', '.')));
        p.MarkerFaceColor = p.Color;
    end
    legend;
    hold off;
    saveas(gcf, sprintf('figures_layout/gvw_rmse_v%s.png', num2str(selected_speed)));
end
