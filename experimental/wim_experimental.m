clc;
close all;
clear;

%% Configurations

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

peak_left = zeros(1, length(pulse_start_l));
peak_right = zeros(1, length(pulse_start_r));

area_left = zeros(1, length(pulse_start_l));
area_right = zeros(1, length(pulse_start_r));

for i = 1:length(pulse_start_l)
    pulse_data = data(pulse_start_l(i):pulse_end_l(i), 1);
    
    % Peak value
    peak_left(i) = max(pulse_data);
    
    % Area under the curve
    area_left(i) = trapz(pulse_data);
    
    % Re-sampling of area (TODO)
end

for i = 1:length(pulse_start_r)
    pulse_data = data(pulse_start_r(i):pulse_end_r(i), 2);
    
    % Peak value
    peak_right(i) = max(pulse_data); % Find the peak value
    
    % Area under the curve
    area_right(i) = trapz(pulse_data); % Calculate the area under the curve
    
    % Re-sampling of area (TODO)
end

%% Display Results
STR = ['Left axle peak: ', num2str(peak_left)];
disp(STR);
STR = ['Right axle peak: ', num2str(peak_right)];
disp(STR);

STR = ['Left axle area: ', num2str(area_left)];
disp(STR);
STR = ['Right axle area: ', num2str(area_right)];
disp(STR);

csvName = ['outputs/results.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'axle,peak_left,peak_right,area_left,area_right\n');
else
    csvFile = fopen(csvName, 'a');
end

for i = 1:3
    fprintf(csvFile, '%d,%.3f,%.3f,%.3f,%.3f\n', i, peak_left(i), peak_right(i), area_left(i), area_right(i));
end