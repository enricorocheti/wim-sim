clc;
close all;
clear;

%% Configurations

% Static configurations
vehicle_axles = 3;
vehicle_passes = 3;
sensor_width = 0.07; % 7 cm

for Fs = [12, 20, 25, 30, 35]
for speed = [5, 10, 15, 20, 30, 40]
for k = 1:vehicle_passes

if Fs == 12 && speed == 5 && k == 1
    continue; % bad csv data
end

if Fs == 12 && speed == 10 && k == 2
    continue; % bad csv data
end

if Fs == 25 && speed == 5 && k == 2
    continue; % bad csv data
end

if Fs == 35 && ((speed == 30 && k == 2) || (speed == 15 && k == 3) || (speed == 15 && k == 2) || (speed == 10 && k == 3) || (speed == 10 && k == 2) || (speed == 5))
    continue; % bad csv data
end

% Load data from CSV
file = ['csv/',num2str(Fs),'khz/',num2str(speed),'km_',num2str(Fs),'k (',num2str(k),').csv'];
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

% figure(2);
% plot(data(pulse_start_l(2):pulse_end_l(2),1));

%% Instantaneous Axle Weight Algorithms

T_resample = sensor_width/(speed/3.6);
samples_per_slice = round((Fs * 1000) * T_resample);

result_peak = zeros(1, vehicle_axles);
result_area = zeros(1, vehicle_axles);
result_resample = zeros(1, vehicle_axles);

for i = 1:vehicle_axles
    pulse_data_left = data(pulse_start_l(i):pulse_end_l(i), 1);
    pulse_data_right = data(pulse_start_r(i):pulse_end_r(i), 2);
    
    % Peak value
    result_peak(i) = max(pulse_data_left) + max(pulse_data_right);
    
    % Area under the curve
    result_area(i) = trapz(pulse_data_left) + trapz(pulse_data_right);
    
    % Re-sampling of area
    % left axle
    num_slices_left = floor(length(pulse_data_left) / samples_per_slice);
    for j = 1:num_slices_left
        start_idx = (j - 1) * samples_per_slice + 1;
        end_idx = min(start_idx + samples_per_slice - 1, length(pulse_data_left));
        result_resample(i) = result_resample(i) + sum(pulse_data_left(start_idx:end_idx));
    end

    % right axle
    num_slices_right = floor(length(pulse_data_right) / samples_per_slice);
    for j = 1:num_slices_right
        start_idx = (j - 1) * samples_per_slice + 1;
        end_idx = min(start_idx + samples_per_slice - 1, length(pulse_data_right));
        result_resample(i) = result_resample(i) + sum(pulse_data_right(start_idx:end_idx));
    end
    
    % Footprint reconstruction (TODO)
    
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

A1_A2_resample = result_resample(1)/result_resample(2);
A1_A3_resample = result_resample(1)/result_resample(3);
A2_A3_resample = result_resample(2)/result_resample(3);

err_A1_A2_peak = (A1_A2_peak - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_peak = (A1_A3_peak - A1_A3_static) * 100 / A1_A3_static;
err_A2_A3_peak = (A2_A3_peak - A2_A3_static) * 100 / A2_A3_static;

err_A1_A2_area = (A1_A2_area - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_area = (A1_A3_area - A1_A3_static) * 100 / A1_A3_static;
err_A2_A3_area = (A2_A3_area - A2_A3_static) * 100 / A2_A3_static;

err_A1_A2_resample = (A1_A2_resample - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_resample = (A1_A3_resample - A1_A3_static) * 100 / A1_A3_static;
err_A2_A3_resample = (A2_A3_resample - A2_A3_static) * 100 / A2_A3_static;

%% Display Results
STR = ['Peak result: ', num2str(result_peak)];
disp(STR);
STR = ['Area result: ', num2str(result_area)];
disp(STR);

%csvName = ['outputs/results',num2str(speed),'.csv'];
csvName = ['outputs/results.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'speed,fs,file_idx,');
    fprintf(csvFile, 'A1/A2_peak,A1/A3_peak,A2/A3_peak,');
    fprintf(csvFile, 'A1/A2_area,A1/A3_area,A2/A3_area,');
    fprintf(csvFile, 'A1/A2_resample,A1/A3_resample,A2/A3_resample\n');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '%d,%dk,(%d),', speed, Fs, k);
fprintf(csvFile, '%.3f,%.3f,%.3f,', err_A1_A2_peak, err_A1_A3_peak, err_A2_A3_peak);
fprintf(csvFile, '%.3f,%.3f,%.3f,', err_A1_A2_area, err_A1_A3_area, err_A2_A3_area);
fprintf(csvFile, '%.3f,%.3f,%.3f\n', err_A1_A2_resample, err_A1_A3_resample, err_A2_A3_resample);
fclose(csvFile);

end
end
end