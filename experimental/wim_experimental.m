clc;
close all;
clear;

%% Configurations

% Static configurations
vehicle_axles = 3;
vehicle_passes = 3;

for Fs = [12, 20, 25, 30]
for k = 1:vehicle_passes

% Load data from CSV
%file = 'csv/25khz/40km_25k (1).csv';
%file = ['csv/',num2str(i),'khz/',num2str(j),'kmh_',num2str(i),'k (',num2str(k),').csv'];
file = ['csv/',num2str(Fs),'khz/40km_',num2str(Fs),'k (',num2str(k),').csv'];
%disp(file);
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

% % Left axle
% figure;
% subplot(2,1,1);
% plot(data(:,1));
% title('Normalized Left Axle');
% hold on;
% plot([1 length(data(:,1))], [threshold_l threshold_l], 'r--');
% hold off;
% 
% % Right axle
% subplot(2,1,2);
% plot(data(:,2));
% title('Normalized Right Axle');
% hold on;
% plot([1 length(data(:,2))], [threshold_r threshold_r], 'r--');
% hold off;

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

%% Errors and statistics
A1_A2_static = 6000/8200;
A1_A3_static = 6000/8200;
A2_A3_static = 8200/8200;

A1_A2_peak = result_peak(1)/result_peak(2);
A1_A3_peak = result_peak(1)/result_peak(3);
A2_A3_peak = result_peak(2)/result_peak(3);

A1_A2_area = result_area(1)/result_area(2);
A1_A3_area = result_area(1)/result_area(3);
A2_A3_area = result_area(2)/result_area(3);

err_A1_A2_peak = (A1_A2_peak - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_peak = (A1_A3_peak - A1_A3_static) * 100 / A1_A3_static;
err_A2_A3_peak = (A2_A3_peak - A2_A3_static) * 100 / A2_A3_static;

err_A1_A2_area = (A1_A2_area - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_area = (A1_A3_area - A1_A3_static) * 100 / A1_A3_static;
err_A2_A3_area = (A2_A3_area - A2_A3_static) * 100 / A2_A3_static;

%% Display Results
STR = ['Peak result: ', num2str(result_peak)];
disp(STR);
STR = ['Area result: ', num2str(result_area)];
disp(STR);

csvName = ['outputs/results.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'speed,fs,file_idx,axle,A1/A2_peak,A1/A3_peak,A2/A3_peak,A1/A2_area,A1/A3_area,A2/A3_area\n');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '40,%dk,(%d),%d,', Fs, k, i);
fprintf(csvFile, '%.3f,%.3f,%.3f,', err_A1_A2_peak, err_A1_A3_peak, err_A2_A3_peak);
fprintf(csvFile, '%.3f,%.3f,%.3f\n', err_A1_A2_area, err_A1_A3_area, err_A2_A3_area);
fclose(csvFile);

end
end