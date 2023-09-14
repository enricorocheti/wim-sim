clc
close all
clear

%% Vehicle parameters
v_speed = 100;               % vehicle speed (km/h)    
v_speed = v_speed/3.6;
%v_acc = 0;                 % FUTURE: vehicle acceleration (m/s^2)
%v_ax_qtty = 2;             % FUTURE: number of axles on vehicle
%v_ax_dist = 1.8;           % FUTURE: distance between those axles

%% System parameters and layout
s_length = 0.07;            % sensor length (meters)
s_qtty = 5;                 % number of sensors
s_dist = 2.5;               % distance between sensors (meters)
s_time = s_length/v_speed;  % time that the sensor is measuring the weight (meters)
s_pos = zeros(1,s_qtty);    % sensor position (meters)
for i = 1:s_qtty
    s_pos(1,i) = i * s_dist;
end

%% Weight signal parameters
w_time_end = s_pos(end)/v_speed;    % time to the vehicle travel through all sensors
w_time_end = w_time_end + 0.1;      % adding some seconds as a margin

w_time_res = 0.000001;                      	% time resolution (seconds)
t = (0:w_time_res:w_time_end-w_time_res);       % time vector

w_f1 = 2;               % first dynamic load frequency (Hz)
w_f1_amp = 15;          % first dynamic load amplitude (Kg)     TODO: real amplitude
w_f2 = 10;              % second dynamic load frequency (Hz)
w_f2_amp = 10;          % second dynamic load amplitude (Kg)    TODO: real amplitude
w_static_load = 100;  	% Axle weight                           TODO: real weight

w_signal = w_static_load + w_f1_amp*sin(2*pi*w_f1*t) + w_f2_amp*sin(2*pi*w_f2*t);

%% ADC sampling
sample_freq = 100e3;                	% ADC sample rate
sample_qtty = s_time/(1/sample_freq);   % Number of samples per sensor

%% Simulations
s_t_ini = zeros(1,s_qtty);      % sensor init time
for i = 1:s_qtty
    s_t_ini(1,i) = s_pos(i)/v_speed;
end

s_idx_ini = zeros(1,s_qtty);    % sensor init index
for i = 1:s_qtty
    s_idx_ini(1,i) = find( abs(t-s_t_ini(i)) < (w_time_res*0.99) );
end

t_size = size(t,2);
s_w_signal = zeros(1,t_size);   % sensor signal

for i = 1:s_qtty
    for j = s_idx_ini(1,i):s_idx_ini(1,i)+sample_qtty
        s_w_signal(j) = w_signal(j);
    end
end

%% Outputs
figure(1);
plot(t,w_signal,'r');
hold on
plot(t,s_w_signal,'b');
xlabel('t [seconds]');
ylabel('Load [Kg]');
line ([t(1) t(end)],[w_static_load w_static_load],'linestyle', '--','color', 'g');
legend('Dynamic axle load','Sensor measures','Static axle load');
title("Vehicle speed = " + v_speed*3.6 + " km/h, sensor distance = " + s_dist + " meters");

s_w_static_load = findpeaks(s_w_signal);
s_w_mean = mean(s_w_static_load);
error = (s_w_mean - w_static_load)*100/w_static_load;