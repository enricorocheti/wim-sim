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
% all_data = [s1_l; s1_r; s2_l; s2_r; s3_l; s3_r];
% common_y_limits = [0, max(all_data)];
% common_x_margin = 10000;
% 
% figure('Name', 'Axle Readings', 'NumberTitle', 'off', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% 
% % S1 Left and Right
% subplot(3, 2, 1);
% plot(s1_l);
% title(['S1 - Normalized Left Axle: ', file]);
% hold on;
% plot([1 length(s1_l)], [threshold_s1_l threshold_s1_l], 'r--');
% hold off;
% ylim(common_y_limits);
% xlim([pulse_ini_s1_l(1)-common_x_margin, pulse_end_s1_l(end)+common_x_margin]);
% 
% subplot(3, 2, 2);
% plot(s1_r);
% title('S1 - Normalized Right Axle');
% hold on;
% plot([1 length(s1_r)], [threshold_s1_r threshold_s1_r], 'r--');
% hold off;
% ylim(common_y_limits);
% xlim([pulse_ini_s1_r(1)-common_x_margin, pulse_end_s1_r(end)+common_x_margin]);
% 
% % S2 Left and Right
% subplot(3, 2, 3);
% plot(s2_l);
% title('S2 - Normalized Left Axle');
% hold on;
% plot([1 length(s2_l)], [threshold_s2_l threshold_s2_l], 'r--');
% hold off;
% ylim(common_y_limits);
% xlim([pulse_ini_s2_l(1)-common_x_margin, pulse_end_s2_l(end)+common_x_margin]);
% 
% subplot(3, 2, 4);
% plot(s2_r);
% title('S2 - Normalized Right Axle');
% hold on;
% plot([1 length(s2_r)], [threshold_s2_r threshold_s2_r], 'r--');
% hold off;
% ylim(common_y_limits);
% xlim([pulse_ini_s2_r(1)-common_x_margin, pulse_end_s2_r(end)+common_x_margin]);
% 
% % S3 Left and Right
% subplot(3, 2, 5);
% plot(s3_l);
% title('S3 - Normalized Left Axle');
% hold on;
% plot([1 length(s3_l)], [threshold_s3_l threshold_s3_l], 'r--');
% hold off;
% ylim(common_y_limits);
% xlim([pulse_ini_s3_l(1)-common_x_margin, pulse_end_s3_l(end)+common_x_margin]);
% 
% subplot(3, 2, 6);
% plot(s3_r);
% title('S3 - Normalized Right Axle');
% hold on;
% plot([1 length(s3_r)], [threshold_s3_r threshold_s3_r], 'r--');
% hold off;
% ylim(common_y_limits);
% xlim([pulse_ini_s3_r(1)-common_x_margin, pulse_end_s3_r(end)+common_x_margin]);

%% Instantaneous Axle Weight Algorithms

T_resample = sensor_width/(speed/3.6);
samples_per_slice = round(Fs * T_resample);

result_peak      = zeros(3, vehicle_axles);
result_area      = zeros(3, vehicle_axles);
result_resample  = zeros(3, vehicle_axles);
result_footprint = zeros(3, vehicle_axles);

mean_peak     = zeros(1, vehicle_axles);
mean_area     = zeros(1, vehicle_axles);
mean_resample = zeros(1, vehicle_axles);
mean_freconst = zeros(1, vehicle_axles);

tire_length_L = zeros(1, vehicle_axles);
tire_length_R = zeros(1, vehicle_axles);
footprint_L   = zeros(1, vehicle_axles);
footprint_R   = zeros(1, vehicle_axles);

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
    % Footprint reconstruction
    tire_length_L(i) = tire_length_L(i) + ((speed/3.6)*(length(pulse_data_left)/Fs));
    tire_length_R(i) = tire_length_R(i) + ((speed/3.6)*(length(pulse_data_right)/Fs));
    
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
    % Footprint reconstruction
    tire_length_L(i) = tire_length_L(i) + ((speed/3.6)*(length(pulse_data_left)/Fs));
    tire_length_R(i) = tire_length_R(i) + ((speed/3.6)*(length(pulse_data_right)/Fs));
    
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
    % Footprint reconstruction
    tire_length_L(i) = tire_length_L(i) + ((speed/3.6)*(length(pulse_data_left)/Fs));
    tire_length_R(i) = tire_length_R(i) + ((speed/3.6)*(length(pulse_data_right)/Fs));
    
    %% Mean S1, S2, S3
    mean_peak(i)     = mean(result_peak(:,i));
    mean_area(i)     = mean(result_area(:,i)) * speed;  % speed compensation
    mean_resample(i) = mean(result_resample(:,i));
    
    % average tire length
    tire_length_L(i) = tire_length_L(i)/3;
    tire_length_R(i) = tire_length_R(i)/3;
    
    if i == 1
        % single tire width is 300 mm
        footprint_L(i) = tire_length_L(i) * 0.3;
        footprint_R(i) = tire_length_R(i) * 0.3;
    else
        % double tire width is 600 mm
        footprint_L(i) = tire_length_L(i) * 0.6;
        footprint_R(i) = tire_length_R(i) * 0.6;
    end
    
    % test this with area and with peak
    mean_freconst(i) = (footprint_L(i) + footprint_R(i)) * mean_area(i); 
    
end

%% Errors, statistics and outputs
GVW_peak     = mean_peak(1) + mean_peak(2) + mean_peak(3) + mean_peak(4) + mean_peak(5) + mean_peak(6);
GVW_area     = mean_area(1) + mean_area(2) + mean_area(3) + mean_area(4) + mean_area(5) + mean_area(6);
GVW_resample = mean_resample(1) + mean_resample(2) + mean_resample(3) + mean_resample(4) + mean_resample(5) + mean_resample(6);
GVW_freconst = mean_freconst(1) + mean_freconst(2) + mean_freconst(3) + mean_freconst(4) + mean_freconst(5) + mean_freconst(6);

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

A1_A2_freconst = mean_freconst(1)/mean_freconst(2);
A1_A3_freconst = mean_freconst(1)/mean_freconst(3);
A1_A4_freconst = mean_freconst(1)/mean_freconst(4);
A1_A5_freconst = mean_freconst(1)/mean_freconst(5);
A1_A6_freconst = mean_freconst(1)/mean_freconst(6);

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

err_A1_A2_freconst = (A1_A2_freconst - A1_A2_static) * 100 / A1_A2_static;
err_A1_A3_freconst = (A1_A3_freconst - A1_A3_static) * 100 / A1_A3_static;
err_A1_A4_freconst = (A1_A4_freconst - A1_A4_static) * 100 / A1_A4_static;
err_A1_A5_freconst = (A1_A5_freconst - A1_A5_static) * 100 / A1_A5_static;
err_A1_A6_freconst = (A1_A6_freconst - A1_A6_static) * 100 / A1_A6_static;

csvName = ['outputs/results_accuracy.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'speed,file_idx,');
    fprintf(csvFile, 'A12_peak,A13_peak,A14_peak,A15_peak,A16_peak,');
    fprintf(csvFile, 'A12_area,A13_area,A14_area,A15_area,A16_area,');
    fprintf(csvFile, 'A12_resample,A13_resample,A14_resample,A15_resample,A16_resample');
    fprintf(csvFile, 'A12_freconst,A13_freconst,A14_freconst,A15_freconst,A16_freconst\n');
    csvFile = fopen(csvName, 'a');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '%d,(%d),', speed, k);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,',	err_A1_A2_peak,     err_A1_A3_peak,     err_A1_A4_peak,     err_A1_A5_peak,     err_A1_A6_peak);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,',	err_A1_A2_area,     err_A1_A3_area,     err_A1_A4_area,     err_A1_A5_area,     err_A1_A6_area);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f',    err_A1_A2_resample, err_A1_A3_resample, err_A1_A4_resample, err_A1_A5_resample, err_A1_A6_resample);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f\n',  err_A1_A2_freconst, err_A1_A3_freconst, err_A1_A4_freconst, err_A1_A5_freconst, err_A1_A6_freconst);
fclose(csvFile);

csvName = ['outputs/results_precision.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'speed,file_idx,');
    fprintf(csvFile, 'A1_peak,A2_peak,A3_peak,A4_peak,A5_peak,A6_peak,');
    fprintf(csvFile, 'A1_area,A2_area,A3_area,A4_area,A5_area,A6_area,');
    fprintf(csvFile, 'A1_resample,A2_resample,A3_resample,A4_resample,A5_resample,A6_resample,');
    fprintf(csvFile, 'A1_freconst,A2_freconst,A3_freconst,A4_freconst,A5_freconst,A6_freconst,');
    fprintf(csvFile, 'GVW_peak,GVW_area,GVW_resample,GVW_freconst\n');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '%d,(%d),', speed, k);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,',  mean_peak(1),       mean_peak(2),       mean_peak(3),       mean_peak(4),       mean_peak(5),       mean_peak(6));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,',  mean_area(1),       mean_area(2),       mean_area(3),       mean_area(4),       mean_area(5),       mean_area(6));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,',  mean_resample(1),   mean_resample(2),   mean_resample(3),   mean_resample(4),   mean_resample(5),   mean_resample(6));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,',  mean_freconst(1),   mean_freconst(2),   mean_freconst(3),   mean_freconst(4),   mean_freconst(5),   mean_freconst(6));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f\n',           GVW_peak,           GVW_area,           GVW_resample,       GVW_freconst);
fclose(csvFile);

end
end

function plotSignalComponents(pulse_data, title_str)
    % Initialize the best odd component (with the highest symmetry)
    best_odd_component = inf; % Start with a large value
    best_shift = 0; % Initialize the best shift to zero
    min_asymmetry = inf; % Start with a large value for the asymmetry measure
    
    % Iterate over possible shifts to find the best symmetry
    for shift = -100:100 % You can adjust this range based on your signal properties
        % Shift the pulse data
        shifted_data = circshift(pulse_data, shift);
        
        % Decompose the signal into even and odd components
        even_component = 0.5 * (shifted_data + flipud(shifted_data));
        odd_component = 0.5 * (shifted_data - flipud(shifted_data));
        
        % Calculate the asymmetry as the sum of absolute values of the odd component
        asymmetry = sum(abs(odd_component));
        
        % Check if this is the minimum asymmetry found so far
        if asymmetry < min_asymmetry
            min_asymmetry = asymmetry;
            best_odd_component = odd_component;
            best_even_component = even_component;
            best_shift = shift;
        end
    end
    
    % Shift the original pulse_data by the best shift found
    pulse_data = circshift(pulse_data, best_shift);
    
    % Decompose the signal into even and odd components with the best shift
    even_component = best_even_component;
    odd_component = best_odd_component;
    
    % Find the maximum and minimum points of the odd component
    [max_value, max_index] = max(odd_component);
    [min_value, min_index] = min(odd_component);
    
    % Create a figure
    figure;
    
    % Plotting the original signal
    subplot(3,1,1); % three rows, first plot
    plot(pulse_data, 'k-'); % Original signal in black
    title([title_str ' - Original Signal (Shifted by ' num2str(best_shift) ')']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    grid on;
    
    % Plotting the even component
    subplot(3,1,2); % three rows, second plot
    plot(even_component, 'b-'); % Even component in blue
    title([title_str ' - Even Component']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    grid on;
    
    % Plotting the odd component
    subplot(3,1,3); % three rows, third plot
    plot(odd_component, 'r-'); % Odd component in red
    hold on; % Hold the plot for additional plotting
    plot(max_index, max_value, 'g*', 'MarkerSize', 10); % Green star for max
    plot(min_index, min_value, 'mo', 'MarkerSize', 10); % Magenta circle for min
    hold off; % Release the plot hold
    title([title_str ' - Odd Component (Minimized Asymmetry)']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    grid on;
end

