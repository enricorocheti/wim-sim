% Parameters
sampling_rate = 1000;   % Sampling rate (samples per second)
duration = 18;          % Duration of the signal in seconds
num_axles = 3;          % Number of axles

% Time vector
t = linspace(0, duration, duration * sampling_rate);

% Simulate signal
signal = zeros(1, length(t));

% Axle information: position and weight
axle_info = [
    3, 0.6;     % Position and weight for first axle
    12, 0.8;    % Position and weight for second axle
    14, 0.85    % Position and weight for third axle
];

% Simulate axle impacts using Gaussian peaks with different weights
for i = 1:num_axles
    position = axle_info(i, 1);
    weight = axle_info(i, 2);
    
    impact = weight * exp(-((t - position).^2) / (2 * 0.2^2)); % Gaussian peak with weight
    signal = signal + impact;
end

% Plot the signal
figure(1);
plot(t, signal,'LineWidth',2);
title('Signal from Weighing Sensor with Axle Impacts');
xlabel('Time (s)');
ylabel('Signal Amplitude');
ylim([-0.1 1]);
