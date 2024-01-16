%clc
close all
clear

%% Vehicle parameters
v_speed = 60;               % vehicle speed (km/h)    
v_speed = v_speed/3.6;
%v_acc = 0;                 % FUTURE: vehicle acceleration (m/s^2)
%v_ax_qtty = 2;             % FUTURE: number of axles on vehicle
%v_ax_dist = 1.8;           % FUTURE: distance between those axles

%% System parameters and layout
s_length = 0.07;            % sensor length (meters)
s_qtty = 16;                 % number of sensors
s_dist = 2;               % distance between sensors (meters)
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
w_f1_amp = 0.8;          % first dynamic load amplitude (Kg)     TODO: real amplitude
w_f2 = 8;              % second dynamic load frequency (Hz)
w_f2_amp = 0.4;          % second dynamic load amplitude (Kg)    TODO: real amplitude
w_static_load = 5;  	% Axle weight                           TODO: real weight

%w_signal = w_static_load + w_f1_amp*sin(2*pi*w_f1*t + pi/2) + w_f2_amp*sin(2*pi*w_f2*t + pi/2);
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
% [left bottom width height]
figureSize = [100, 100, 700, 500];
figure('Position',figureSize);
plot(t,w_signal,'color',[0.8500 0.3250 0.0980],'LineWidth',1.7);
hold on
line ([t(1) t(end)],[w_static_load w_static_load],'linestyle', '--','color', [0 0.4470 0.7410],'LineWidth',1);
for i = 1:s_qtty
    hold on 
    scatter(t(s_idx_ini(i)),s_w_signal(s_idx_ini(i)),'x','MarkerEdgeColor',[0.3 0.3 0.3],'LineWidth',1.7);
end
hold on
plot(t,s_w_signal,'color',[0.3 0.3 0.3]);
xlabel('Time [s]');
ylabel('Load [t]');
legend('Dynamic axle load','Static axle load','Sensor measures');
xlim([0 max(t)])
%title("Vehicle speed = " + v_speed*3.6 + " km/h, sensor distance = " + s_dist + " meters");

s_w_static_load = findpeaks(s_w_signal);
s_w_mean = mean(s_w_static_load);
error = (s_w_mean - w_static_load)*100/w_static_load;


STR = [' Estimated load = ',num2str(s_w_mean)];
disp(STR)
STR = [' Relative error = ',num2str(error)];
disp(STR)
STR = [' N = ',num2str(s_qtty),' sensors'];
disp(STR)
disp(' ')
%STR = [' V = ',num2str(v_speed*3.6),' km/h'];

close all
%disp(STR)