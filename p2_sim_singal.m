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
s_qtty = 4;                 % number of sensors
s_dist = 2.5;               % distance between sensors (meters)
s_pos = zeros(1,s_qtty);    % sensor position (meters)
for i = 1:s_qtty
    s_pos(1,i) = i * s_dist;
end

%% Weight signal parameters
w_time_end = s_pos(end)/v_speed;    % time to the vehicle travel through all sensors
w_time_end = w_time_end + 0.1;      % adding some seconds as a margin

w_time_res = 0.00001;                      	% time resolution (seconds)
t = (0:w_time_res:w_time_end-w_time_res);       % time vector

w_f1 = 2;               % first dynamic load frequency (Hz)
w_f1_amp = 15;          % first dynamic load amplitude (Kg)     TODO: real amplitude
w_f2 = 10;              % second dynamic load frequency (Hz)
w_f2_amp = 10;          % second dynamic load amplitude (Kg)    TODO: real amplitude
w_static_load = 100;  	% Axle weight                           TODO: real weight

w_signal = w_static_load + w_f1_amp*sin(2*pi*w_f1*t) + w_f2_amp*sin(2*pi*w_f2*t);

%% Simulations

% sensor measure time
s_time = zeros(1,s_qtty);      
for i = 1:s_qtty
    s_time(1,i) = s_pos(i)/v_speed;
end

% sensor measure index
s_idx_ini = zeros(1,s_qtty);
for i = 1:s_qtty
    s_idx_ini(1,i) = find( abs(t-s_time(i)) < (w_time_res*0.99) );
end

% sensor signal
t_size = size(t,2);
s_w_signal = zeros(1,t_size);

for i = 1:s_qtty
    for j = s_idx_ini(1,i):s_idx_ini(1,i)+1
        s_w_signal(j) = w_signal(j);
    end
end
s_w_measures = findpeaks(s_w_signal);

%% Outputs
figure(1);
plot(t,w_signal,'r');
hold on
plot(t,s_w_signal,'b');
for i = 1:length(s_w_measures)
    hold on
    plot(s_time(i),s_w_measures(i),'b*')
end
xlabel('t [seconds]');
ylabel('Load [Kg]');
line ([t(1) t(end)],[w_static_load w_static_load],'linestyle', '--','color', 'g');
%legend('Dynamic axle load','Sensor measures','Static axle load');
title("Vehicle speed = " + v_speed*3.6 + " km/h, sensor distance = " + s_dist + " meters");

s_w_mean = mean(s_w_measures);
error = (s_w_mean - w_static_load)*100/w_static_load;

%% Show results
STR = [num2str(s_qtty),' sensors, speed = ',num2str(v_speed*3.6),' km/h'];
disp(STR)
STR = ['f1 = ',num2str(w_f1),' Hz, f2 = ',num2str(w_f2),' Hz'];
disp(STR)
STR = ['=================================='];
disp(STR)
STR = ['static load = ',num2str(w_static_load),' Kg'];
disp(STR)
STR = ['estim. load = ',num2str(s_w_mean),' Kg'];
disp(STR)
STR = ['error = ',num2str(error),'%'];
disp(STR)