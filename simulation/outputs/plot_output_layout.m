clc
close all
clear

%% Configuration
sensor_dist = 0.5;
sensor_counts = [2, 4, 6, 8, 10, 12, 14, 16];
markers = {'-o', '-s', '-^', '-d', '-p', '-h', '-x', '-+'};

axle_data = cell(length(sensor_counts), 1);
gvw_data = cell(length(sensor_counts), 1);

% Load data from CSV files
for i = 1:length(sensor_counts)
    sensors = sensor_counts(i);
    %axle_data{i} = readtable(sprintf('axle_output_s%d_d%s.csv', sensors, num2str(sensor_dist)));
    %gvw_data{i} = readtable(sprintf('gvw_output_s%d_d%s.csv', sensors, num2str(sensor_dist)));
    axle_data{i} = readtable(sprintf('axle_output_s%d_Delta1.csv', sensors));
    gvw_data{i} = readtable(sprintf('gvw_output_s%d_Delta1.csv', sensors));
end

%% Axle maximum relative error
figure;
hold on;
%title(sprintf('Axle maximum relative error, d = %sm',num2str(sensor_dist)));
title(sprintf('Axle maximum relative error, d = %s','Delta1'));
xlabel('Number of Sensors');
ylabel('Er (%)');
grid on;

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
%saveas(gcf, sprintf('axle_max_er_d%s.png',num2str(sensor_dist)));
saveas(gcf, 'figures_layouts/axle_max_er_Delta1.png');

%% GVW maximum relative error
figure;
hold on;
%title(sprintf('GVW maximum relative error, d = %sm',num2str(sensor_dist)));
title(sprintf('GVW maximum relative error, d = %s','Delta1'));
xlabel('Number of Sensors');
ylabel('Er (%)');
grid on;

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
%saveas(gcf, sprintf('gvw_max_er_d%s.png',num2str(sensor_dist)));
saveas(gcf, 'figures_layouts/gvw_max_er_Delta1.png');

%% Axle absolute mean error
figure;
hold on;
%title(sprintf('Axle absolute mean error, d = %sm',num2str(sensor_dist)));
title(sprintf('Axle absolute mean error, d = %s','Delta1'));
xlabel('Number of Sensors');
ylabel('Mean error (%)');
grid on;

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    mean_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = axle_data{i};
        mean_mv_values(i) = df.mean_mv(df.speed == speed);
    end
    p = plot(sensor_counts, mean_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
%saveas(gcf, sprintf('axle_abs_mean_d%s.png',num2str(sensor_dist)));
saveas(gcf, 'figures_layouts/axle_abs_mean_Delta1.png');

%% GVW absolute mean error
figure;
hold on;
%title(sprintf('GVW absolute mean error, d = %sm',num2str(sensor_dist)));
title(sprintf('GVW absolute mean error, d = %s','Delta1'));
xlabel('Number of Sensors');
ylabel('Mean error (%)');
grid on;

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    mean_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = gvw_data{i};
        mean_mv_values(i) = df.mean_mv(df.speed == speed);
    end
    p = plot(sensor_counts, mean_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
%saveas(gcf, sprintf('gvw_abs_mean_d%s.png',num2str(sensor_dist)));
saveas(gcf, 'figures_layouts/gvw_abs_mean_Delta1.png');

%% Axle standard deviation
figure;
hold on;
%title(sprintf('Axle standard deviation, d = %sm',num2str(sensor_dist)));
title(sprintf('Axle standard deviation, d = %s','Delta1'));
xlabel('Number of Sensors');
ylabel('Standard deviation');
grid on;

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    std_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = axle_data{i};
        std_mv_values(i) = df.std_mv(df.speed == speed);
    end
    p = plot(sensor_counts, std_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
%saveas(gcf, sprintf('axle_std_d%s.png',num2str(sensor_dist)));
saveas(gcf, 'figures_layouts/axle_std_Delta1.png');

%% GVW standard deviation
figure;
hold on;
%title(sprintf('GVW standard deviation, d = %sm',num2str(sensor_dist)));
title(sprintf('GVW standard deviation, d = %s','Delta1'));
xlabel('Number of Sensors');
ylabel('Standard deviation');
grid on;

unique_speeds = unique(axle_data{1}.speed);
for speed_idx = 1:length(unique_speeds)
    speed = unique_speeds(speed_idx);
    std_mv_values = zeros(length(sensor_counts), 1);
    for i = 1:length(sensor_counts)
        df = gvw_data{i};
        std_mv_values(i) = df.std_mv(df.speed == speed);
    end
    p = plot(sensor_counts, std_mv_values, markers{speed_idx}, 'DisplayName', sprintf('Speed %d km/h', speed));
    p.MarkerFaceColor = p.Color;
end
legend;
hold off;
%saveas(gcf, sprintf('gvw_std_d%s.png',num2str(sensor_dist)));
saveas(gcf, 'figures_layouts/gvw_std_Delta1.png');
