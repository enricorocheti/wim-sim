clc
close all
clear

%% Configurations

% Load data from CSV
data = csvread('20km_20k_1.csv');

% Normalize the data
data(:,1) = data(:,1) - min(data(:,1));
data(:,2) = data(:,2) - min(data(:,2));

%% Peak value

% Left axle
figure;
subplot(2,1,1);
plot(data(:,1));
title('Normalized Left Axle');
hold on;

[pks_L, locs_L] = findpeaks(data(:,1), 'MinPeakHeight', 3*std(data(:,1)));
plot(locs_L, pks_L, 'vr');

% Right axle
subplot(2,1,2);
plot(data(:,2));
title('Normalized Right Axle');
hold on
[pks_R, locs_R] = findpeaks(data(:,2), 'MinPeakHeight', 3*std(data(:,2)));
plot(locs_R, pks_R, 'vg');

%% Area under the signal
area_L = zeros(length(locs_L),1);
area_R = zeros(length(locs_R),1);
window = 50;

for i = 1:length(locs_L)
    start_idx = max(1, locs_L(i) - window);
    end_idx = min(length(data(:,1)), locs_L(i) + window);
    area_L(i) = trapz(data(start_idx:end_idx, 1));
end

for i = 1:length(locs_R)
    start_idx = max(1, locs_R(i) - window);
    end_idx = min(length(data(:,2)), locs_R(i) + window);
    area_R(i) = trapz(data(start_idx:end_idx, 2));
end

disp('Areas under peaks for Left Axle:');
disp(sum(area_L));
disp('Areas under peaks for Right Axle:');
disp(sum(area_R));

%% Re-sampling of area
% TODO

%% Footprint reconstruction
% TODO

