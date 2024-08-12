clc;
close all;
clear;

%% Configurations

Fs = 25400;
vehicle_axles = 6;
vehicle_passes = 1;
sensor_width = 0.07;    % 7 cm
load_sens = 3;          % S1, S2, S3
threshold = 0.05;       % 5%

% option 1: loaded truck files
% option 2: empty truck files
option = 1;

if option == 1
    %speeds = [5, 10, 15, 20, 25, 30, 35];
    speeds = [30];
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

% Remove corrupted data from analysis
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

sample_limit = 20000;

% Remove first few and last samples
s1_l = s1_l(sample_limit:end-sample_limit);
s1_r = s1_r(sample_limit:end-sample_limit);
s2_l = s2_l(sample_limit:end-sample_limit);
s2_r = s2_r(sample_limit:end-sample_limit);
s3_l = s3_l(sample_limit:end-sample_limit);
s3_r = s3_r(sample_limit:end-sample_limit);

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
common_y_limits = [0, max(all_data)+5000];
common_x_margin = 10000;

% Figure size
fig_x = 50;
fig_y = 50;
fig_w = 540;
fig_h = 400;
fontSize = 15;
fontSize2 = 12;


%% S1 left and right
figure('Position', [fig_x, fig_y, fig_w, fig_h]);
plot(s1_l, 'DisplayName', 'Signal');
hold on;
plot([1 length(s1_l)], [threshold_s1_l threshold_s1_l], 'r--', 'DisplayName', 'Threshold');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s1_l(1)-common_x_margin, pulse_end_s1_l(end)+common_x_margin]);
set(gca, 'FontSize', fontSize2);
xlabel('Samples', 'FontSize', fontSize);
ylabel('ADC counts', 'FontSize', fontSize);
saveas(gcf, sprintf('outputs/figures/S1L_v%s_%s.png', num2str(speed), num2str(k)));

figure('Position', [fig_x, fig_y, fig_w, fig_h]);
plot(s1_r, 'DisplayName', 'Signal');
hold on;
plot([1 length(s1_r)], [threshold_s1_r threshold_s1_r], 'r--', 'DisplayName', 'Threshold');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s1_r(1)-common_x_margin, pulse_end_s1_r(end)+common_x_margin]);
set(gca, 'FontSize', fontSize2);
xlabel('Samples', 'FontSize', fontSize);
ylabel('ADC counts', 'FontSize', fontSize);
saveas(gcf, sprintf('outputs/figures/S1R_v%s_%s.png', num2str(speed), num2str(k)));

%% S2 left and right
figure('Position', [fig_x, fig_y, fig_w, fig_h]);
plot(s2_l, 'DisplayName', 'Signal');
hold on;
plot([1 length(s2_l)], [threshold_s2_l threshold_s2_l], 'r--', 'DisplayName', 'Threshold');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s2_l(1)-common_x_margin, pulse_end_s2_l(end)+common_x_margin]);
set(gca, 'FontSize', fontSize2);
xlabel('Samples', 'FontSize', fontSize);
ylabel('ADC counts', 'FontSize', fontSize);
saveas(gcf, sprintf('outputs/figures/S2L_v%s_%s.png', num2str(speed), num2str(k)));

figure('Position', [fig_x, fig_y, fig_w, fig_h]);
plot(s2_r, 'DisplayName', 'Signal');
hold on;
plot([1 length(s2_r)], [threshold_s2_r threshold_s2_r], 'r--', 'DisplayName', 'Threshold');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s2_r(1)-common_x_margin, pulse_end_s2_r(end)+common_x_margin]);
set(gca, 'FontSize', fontSize2);
xlabel('Samples', 'FontSize', fontSize);
ylabel('ADC counts', 'FontSize', fontSize);
saveas(gcf, sprintf('outputs/figures/S2R_v%s_%s.png', num2str(speed), num2str(k)));

%% S3 left and right
figure('Position', [fig_x, fig_y, fig_w, fig_h]);
plot(s3_l, 'DisplayName', 'Signal');
hold on;
plot([1 length(s3_l)], [threshold_s3_l threshold_s3_l], 'r--', 'DisplayName', 'Threshold');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s3_l(1)-common_x_margin, pulse_end_s3_l(end)+common_x_margin]);
set(gca, 'FontSize', fontSize2);
xlabel('Samples', 'FontSize', fontSize);
ylabel('ADC counts', 'FontSize', fontSize);
saveas(gcf, sprintf('outputs/figures/S3L_v%s_%s.png', num2str(speed), num2str(k)));

figure('Position', [fig_x, fig_y, fig_w, fig_h]);
plot(s3_r, 'DisplayName', 'Signal');
hold on;
plot([1 length(s3_r)], [threshold_s3_r threshold_s3_r], 'r--', 'DisplayName', 'Threshold');
hold off;
ylim(common_y_limits);
xlim([pulse_ini_s3_r(1)-common_x_margin, pulse_end_s3_r(end)+common_x_margin]);
set(gca, 'FontSize', fontSize2);
xlabel('Samples', 'FontSize', fontSize);
ylabel('ADC counts', 'FontSize', fontSize);
saveas(gcf, sprintf('outputs/figures/S3R_v%s_%s.png', num2str(speed), num2str(k)));

end
end