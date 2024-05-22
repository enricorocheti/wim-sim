clc;
close all;
clear;

%% Configurations

% Static configurations
vehicle_axles = 3;
vehicle_passes = 3;

% Load data from CSV
file = 'csv/25khz/40km_25k (1).csv';
data = csvread(file);

% Normalize the data
data(:,1) = data(:,1) - min(data(:,1)); % left axle
data(:,2) = data(:,2) - min(data(:,2)); % right axle

% Moving average filter
window_size = 5;
data(:,1) = movmean(data(:,1), window_size);
data(:,2) = movmean(data(:,2), window_size);

% Start and end points of pulses using 2.5% threshold
threshold_l = min(data(:,1)) + 0.025 * (max(data(:,1)) - min(data(:,1)));
pulse_start_l = find(diff(data(:,1) > threshold_l) == 1);
pulse_end_l   = find(diff(data(:,1) > threshold_l) == -1);

threshold_r = min(data(:,2)) + 0.025 * (max(data(:,2)) - min(data(:,2)));
pulse_start_r = find(diff(data(:,2) > threshold_r) == 1);
pulse_end_r   = find(diff(data(:,2) > threshold_r) == -1);

%% Plot

% Left axle
figure;
subplot(2,1,1);
plot(data(:,1));
title('Normalized Left Axle');
hold on;
plot([1 length(data(:,1))], [threshold_l threshold_l], 'r--');
hold off;

% Right axle
subplot(2,1,2);
plot(data(:,2));
title('Normalized Right Axle');
hold on;
plot([1 length(data(:,2))], [threshold_r threshold_r], 'r--');
hold off;

%figure(2);
%plot(data(pulse_start_l(2):pulse_end_l(2),1));

%% Instantaneous Axle Weight Algorithms

result_peak = zeros(1, vehicle_axles);
result_area = zeros(1, vehicle_axles);

for i = 1:vehicle_axles
    pulse_data_left = data(pulse_start_l(i):pulse_end_l(i), 1);
    pulse_data_right = data(pulse_start_r(i):pulse_end_r(i), 2);
    
    % Peak value
    result_peak(i) = max(pulse_data_left) + max(pulse_data_right);
    
    % Area under the curve
    result_area(i) = trapz(pulse_data_left) + trapz(pulse_data_right);
    
    % Re-sampling of area (TODO)
end

%% Display Results
STR = ['Peak result: ', num2str(result_peak)];
disp(STR);
STR = ['Area result: ', num2str(result_area)];
disp(STR);

csvName = ['outputs/results2.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'axle,peak,area\n');
else
    csvFile = fopen(csvName, 'a');
end

for i = 1:vehicle_axles
    fprintf(csvFile, '%d,%.3f,%.3f\n', i, result_peak(i), result_area(i));
end