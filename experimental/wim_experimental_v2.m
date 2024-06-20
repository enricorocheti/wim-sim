clc;
close all;
clear;

%% Configurations

Fs = 25400;
vehicle_axles = 6;
vehicle_passes = 10;
sensor_width = 0.07;    % 7 cm
load_sens = 3;          % S1, S2, S3
threshold = 0.05;       % 5%

% option 1: loaded truck files
% option 2: empty truck files
option = 1;

if option == 1
    speeds = [5, 10, 15, 20, 25, 30, 35];
    A1_static = 5957.5;
    A2_static = 9092.5;
    A3_static = 6182.5;
    A4_static = 6960;
    A5_static = 9157.5;
    A6_static = 7532.5;
    GVW_static = 44460;
else
    speeds = [10, 20, 30, 40, 50, 60, 65];
    A1_static = 5630;
    A2_static = 3767.5;
    A3_static = 2370;
    A4_static = 1827.5;
    A5_static = 2860;
    A6_static = 2237.5;
    GVW_static = 18380;
end

for speed = speeds
for k = 1:vehicle_passes

if option == 1 && (speed == 10 && k == 2)
    continue; % bad csv data
end
if option == 1 && (speed == 15 && (k == 9 || k == 10))
    continue; % bad csv data
end
if option == 2 && (speed == 30 && (k == 2 || k == 8))
    continue; % bad csv data
end
if option == 2 && (speed == 40 && k == 2)
    continue; % bad csv data
end
if option == 2 && (speed == 60 && k == 5)
    continue; % bad csv data
end
if option == 2 && (speed == 65 && k > 5)
    continue; % there is only 5 files in this speed
end

%% Signal Processing

% Load data from CSV
if option == 1
    file = ['csv/6_axle_truck/Vel_',num2str(speed),'km/',num2str(speed),'km (',num2str(k),').csv'];
else
    file = ['csv/6_axle_truck_empty/Vel_',num2str(speed),'km/',num2str(speed),'km (',num2str(k),').csv'];
end
disp(file);
data = readtable(file);

s1_l = data.ch0_flt;    % S1 left
s1_r = data.ch2_flt;    % S1 right
s2_l = data.ch1_flt;    % S2 left
s2_r = data.ch3_flt;    % S2 right
s3_l = data.ch4_flt;    % S3 left
s3_r = data.ch5_flt;	% S3 right

% Remove first few and last samples
s1_l = s1_l(101:end-100);
s1_r = s1_r(101:end-100);
s2_l = s2_l(101:end-100);
s2_r = s2_r(101:end-100);
s3_l = s3_l(101:end-100);
s3_r = s3_r(101:end-100);

% Normalize the data
s1_l = s1_l - min(s1_l);
s1_r = s1_r - min(s1_r);
s2_l = s2_l - min(s2_l);
s2_r = s2_r - min(s2_r);
s3_l = s3_l - min(s3_l);
s3_r = s3_r - min(s3_r);

% Calculate the start and end of each pulse using threshold (in %)
threshold_s1_l = min(s1_l) + threshold * (max(s1_l) - min(s1_l));
threshold_s2_l = min(s2_l) + threshold * (max(s2_l) - min(s2_l));
threshold_s3_l = min(s3_l) + threshold * (max(s3_l) - min(s3_l));
pulse_ini_s1_l = find(diff(s1_l > threshold_s1_l) == 1);
pulse_ini_s2_l = find(diff(s2_l > threshold_s2_l) == 1);
pulse_ini_s3_l = find(diff(s3_l > threshold_s3_l) == 1);
pulse_end_s1_l = find(diff(s1_l > threshold_s1_l) == -1);
pulse_end_s2_l = find(diff(s2_l > threshold_s2_l) == -1);
pulse_end_s3_l = find(diff(s3_l > threshold_s3_l) == -1);

threshold_s1_r = min(s1_r) + threshold * (max(s1_r) - min(s1_r));
threshold_s2_r = min(s2_r) + threshold * (max(s2_r) - min(s2_r));
threshold_s3_r = min(s3_r) + threshold * (max(s3_r) - min(s3_r));
pulse_ini_s1_r = find(diff(s1_r > threshold_s1_r) == 1);
pulse_ini_s2_r = find(diff(s2_r > threshold_s2_r) == 1);
pulse_ini_s3_r = find(diff(s3_r > threshold_s3_r) == 1);
pulse_end_s1_r = find(diff(s1_r > threshold_s1_r) == -1);
pulse_end_s2_r = find(diff(s2_r > threshold_s2_r) == -1);
pulse_end_s3_r = find(diff(s3_r > threshold_s3_r) == -1);

%% Figures
all_data = [s1_l; s1_r; s2_l; s2_r; s3_l; s3_r];
common_y_limits = [0, max(all_data)];
common_x_margin = 10000;

figure('Name', 'Axle Readings', 'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% S1 Left and Right
subplot(3, 2, 1);
plot(s1_l);
title(['S1 - Normalized Left Axle: ', file]);
hold on;
plot([1 length(s1_l)], [threshold_s1_l threshold_s1_l], 'r--');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s1_l(1)-common_x_margin, pulse_end_s1_l(end)+common_x_margin]);

subplot(3, 2, 2);
plot(s1_r);
title('S1 - Normalized Right Axle');
hold on;
plot([1 length(s1_r)], [threshold_s1_r threshold_s1_r], 'r--');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s1_r(1)-common_x_margin, pulse_end_s1_r(end)+common_x_margin]);

% S2 Left and Right
subplot(3, 2, 3);
plot(s2_l);
title('S2 - Normalized Left Axle');
hold on;
plot([1 length(s2_l)], [threshold_s2_l threshold_s2_l], 'r--');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s2_l(1)-common_x_margin, pulse_end_s2_l(end)+common_x_margin]);

subplot(3, 2, 4);
plot(s2_r);
title('S2 - Normalized Right Axle');
hold on;
plot([1 length(s2_r)], [threshold_s2_r threshold_s2_r], 'r--');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s2_r(1)-common_x_margin, pulse_end_s2_r(end)+common_x_margin]);

% S3 Left and Right
subplot(3, 2, 5);
plot(s3_l);
title('S3 - Normalized Left Axle');
hold on;
plot([1 length(s3_l)], [threshold_s3_l threshold_s3_l], 'r--');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s3_l(1)-common_x_margin, pulse_end_s3_l(end)+common_x_margin]);

subplot(3, 2, 6);
plot(s3_r);
title('S3 - Normalized Right Axle');
hold on;
plot([1 length(s3_r)], [threshold_s3_r threshold_s3_r], 'r--');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s3_r(1)-common_x_margin, pulse_end_s3_r(end)+common_x_margin]);

%% Instantaneous Axle Weight Algorithms

T_resample = sensor_width/(speed/3.6);
samples_per_slice = round(Fs * T_resample);

result_peak = zeros(3, vehicle_axles);
result_area = zeros(3, vehicle_axles);
result_resample = zeros(3, vehicle_axles);
result_footprint = zeros(3, vehicle_axles);

mean_peak = zeros(1, vehicle_axles);
mean_area = zeros(1, vehicle_axles);
mean_resample = zeros(1, vehicle_axles);
mean_footprint = zeros(1, vehicle_axles);

for i = 1:vehicle_axles
    %% S1
    pulse_data_left = s1_l(pulse_ini_s1_l(i):pulse_end_s1_l(i));
    pulse_data_right = s1_r(pulse_ini_s1_r(i):pulse_end_s1_r(i));
    
    % Peak
    result_peak(1,i) = max(pulse_data_left) + max(pulse_data_right);
    % Area
    result_area(1,i) = trapz(pulse_data_left) + trapz(pulse_data_right);
    % Re-sampling of area
    num_slices_left = floor(length(pulse_data_left) / samples_per_slice);
    num_slices_right = floor(length(pulse_data_right) / samples_per_slice);
    for j = 1:num_slices_left
        slice_idx = (j - 1) * samples_per_slice + 1;
        result_resample(1,i) = result_resample(1,i) + pulse_data_left(slice_idx);
    end
    for j = 1:num_slices_right
        slice_idx = (j - 1) * samples_per_slice + 1;
        result_resample(1,i) = result_resample(1,i) + pulse_data_right(slice_idx);
    end
    
    %% S2
    pulse_data_left = s2_l(pulse_ini_s2_l(i):pulse_end_s2_l(i));
    pulse_data_right = s2_r(pulse_ini_s2_r(i):pulse_end_s2_r(i));
    
    % Peak
    result_peak(2,i) = max(pulse_data_left) + max(pulse_data_right);
    % Area
    result_area(2,i) = trapz(pulse_data_left) + trapz(pulse_data_right);
    % Re-sampling of area
    num_slices_left = floor(length(pulse_data_left) / samples_per_slice);
    num_slices_right = floor(length(pulse_data_right) / samples_per_slice);
    for j = 1:num_slices_left
        slice_idx = (j - 1) * samples_per_slice + 1;
        result_resample(2,i) = result_resample(2,i) + pulse_data_left(slice_idx);
    end
    for j = 1:num_slices_right
        slice_idx = (j - 1) * samples_per_slice + 1;
        result_resample(2,i) = result_resample(2,i) + pulse_data_right(slice_idx);
    end
    
    %% S3
    pulse_data_left = s3_l(pulse_ini_s3_l(i):pulse_end_s3_l(i));
    pulse_data_right = s3_r(pulse_ini_s3_r(i):pulse_end_s3_r(i));
    
    % Peak
    result_peak(3,i) = max(pulse_data_left) + max(pulse_data_right);
    % Area
    result_area(3,i) = trapz(pulse_data_left) + trapz(pulse_data_right);
    % Re-sampling of area
    num_slices_left = floor(length(pulse_data_left) / samples_per_slice);
    num_slices_right = floor(length(pulse_data_right) / samples_per_slice);
    for j = 1:num_slices_left
        slice_idx = (j - 1) * samples_per_slice + 1;
        result_resample(3,i) = result_resample(3,i) + pulse_data_left(slice_idx);
    end
    for j = 1:num_slices_right
        slice_idx = (j - 1) * samples_per_slice + 1;
        result_resample(3,i) = result_resample(3,i) + pulse_data_right(slice_idx);
    end
    
    %% Mean S1, S2, S3
    mean_peak(i) = mean(result_peak(:,i));
    mean_area(i) = speed * mean(result_area(:,i));   % speed compensation
    mean_resample(i) = mean(result_resample(:,i));
end

%% Errors and statistics
A1_A2_static = A1_static/A2_static;
A1_A3_static = A1_static/A3_static;
A1_A4_static = A1_static/A4_static;
A1_A5_static = A1_static/A5_static;
A1_A6_static = A1_static/A6_static;

A1_A2_peak = mean_peak(1)/mean_peak(2);
A1_A3_peak = mean_peak(1)/mean_peak(3);
A1_A4_peak = mean_peak(1)/mean_peak(4);
A1_A5_peak = mean_peak(1)/mean_peak(5);
A1_A6_peak = mean_peak(1)/mean_peak(6);

A1_A2_area = mean_area(1)/mean_area(2);
A1_A3_area = mean_area(1)/mean_area(3);
A1_A4_area = mean_area(1)/mean_area(4);
A1_A5_area = mean_area(1)/mean_area(5);
A1_A6_area = mean_area(1)/mean_area(6);

A1_A2_resample = mean_resample(1)/mean_resample(2);
A1_A3_resample = mean_resample(1)/mean_resample(3);
A1_A4_resample = mean_resample(1)/mean_resample(4);
A1_A5_resample = mean_resample(1)/mean_resample(5);
A1_A6_resample = mean_resample(1)/mean_resample(6);

err_A1_A2_peak = (A1_A2_peak - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_peak = (A1_A3_peak - A1_A3_static) * 100 / A1_A3_static;
err_A1_A4_peak = (A1_A4_peak - A1_A4_static) * 100 / A1_A4_static;
err_A1_A5_peak = (A1_A5_peak - A1_A5_static) * 100 / A1_A5_static;
err_A1_A6_peak = (A1_A6_peak - A1_A6_static) * 100 / A1_A6_static;

err_A1_A2_area = (A1_A2_area - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_area = (A1_A3_area - A1_A3_static) * 100 / A1_A3_static;
err_A1_A4_area = (A1_A4_area - A1_A4_static) * 100 / A1_A4_static;
err_A1_A5_area = (A1_A5_area - A1_A5_static) * 100 / A1_A5_static;
err_A1_A6_area = (A1_A6_area - A1_A6_static) * 100 / A1_A6_static;

err_A1_A2_resample = (A1_A2_resample - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_resample = (A1_A3_resample - A1_A3_static) * 100 / A1_A3_static;
err_A1_A4_resample = (A1_A4_resample - A1_A4_static) * 100 / A1_A4_static;
err_A1_A5_resample = (A1_A5_resample - A1_A5_static) * 100 / A1_A5_static;
err_A1_A6_resample = (A1_A6_resample - A1_A6_static) * 100 / A1_A6_static;

GVW_peak = mean_peak(1) + mean_peak(2) + mean_peak(3) + mean_peak(4) + mean_peak(5) + mean_peak(6);
GVW_area = mean_area(1) + mean_area(2) + mean_area(3) + mean_area(4) + mean_area(5) + mean_area(6);
GVW_resample = mean_resample(1) + mean_resample(2) + mean_resample(3) + mean_resample(4) + mean_resample(5) + mean_resample(6);

%% Results
csvName = ['outputs/results_v2.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'speed,file_idx,');
%     fprintf(csvFile, 'A1/A2_peak,A1/A3_peak,A1/A4_peak,A1/A5_peak,A1/A6_peak,');
%     fprintf(csvFile, 'A1/A2_area,A1/A3_area,A1/A4_area,A1/A5_area,A1/A6_area,');
%     fprintf(csvFile, 'A1/A2_resample,A1/A3_resample,A1/A4_resample,A1/A5_resample,A1/A6_resample,');
    fprintf(csvFile, 'A1_peak,A2_peak,A3_peak,A4_peak,A5_peak,A6_peak,');
    fprintf(csvFile, 'A1_area,A2_area,A3_area,A4_area,A5_area,A6_area,');
    fprintf(csvFile, 'A1_resample,A2_resample,A3_resample,A4_resample,A5_resample,A6_resample,');
    fprintf(csvFile, 'GVW_peak,GVW_area,GVW_resample\n');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '%d,(%d),', speed, k);
% fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,', err_A1_A2_peak, err_A1_A3_peak, err_A1_A4_peak, err_A1_A5_peak, err_A1_A6_peak);
% fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,', err_A1_A2_area, err_A1_A3_area, err_A1_A4_area, err_A1_A5_area, err_A1_A6_area);
% fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,', err_A1_A2_resample, err_A1_A3_resample, err_A1_A4_resample, err_A1_A5_resample, err_A1_A6_resample);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,',  mean_peak(1),  mean_peak(2),  mean_peak(3),  mean_peak(4),  mean_peak(5),  mean_peak(6));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,', mean_area(1), mean_area(2), mean_area(3), mean_area(4), mean_area(5), mean_area(6));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,', mean_resample(1), mean_resample(2), mean_resample(3), mean_resample(4), mean_resample(5), mean_resample(6));
fprintf(csvFile, '%.3f,%.3f,%.3f\n', GVW_peak, GVW_area, GVW_resample);
fclose(csvFile);

end
end