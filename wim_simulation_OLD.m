clc
close all
clear

%% Vehicle parameters
v_speed = 60;              % vehicle speed (km/h)    
v_speed = v_speed/3.6;

%% Vehicles
v1_axle_qtty = 3;                   % number of axles on vehicle
v1_axle_dist = [4.8, 1.2];          % distance between axles
v1_axle_st_load = [6, 8.2, 8.2];    % axle weight

% TODO: add more vehicles and use a general matrix for that

%% System parameters and layout
s_qtty = 6;                 % number of sensors
s_dist = 2.5;               % distance between sensors (meters)
s_pos = zeros(1,s_qtty);    % sensor position (meters)
for i = 1:s_qtty
    s_pos(1,i) = i * s_dist; % TODO: s_pos starting on zero? now it aint
end

%% Weight signal parameters
w_time_end = s_pos(end)/v_speed;    % time to the vehicle travel through all sensors
w_time_end = w_time_end + 0.1;      % adding time as a margin

w_time_res = 0.001;                           % time resolution (seconds)
t = (0:w_time_res:w_time_end-w_time_res);       % time vector
t_size = size(t,2);                             % time vector size

% f1 and f2 range (Hz)
w_f1_min = 1;               
w_f1_max = 4;
w_f2_min = 8;
w_f2_max = 15;

% TODO: amplitude should to depend on vehicle speed
w_f1_amp = 15;          % first dynamic load amplitude (Kg)    
w_f2_amp = 10;          % second dynamic load amplitude (Kg)

w_static_load = 100;    % TODO: use vX_ax1_st_load instead of this

n = 1000;                % number of simulations

w_signal = zeros(n,t_size);
for i = 1:n
	w_f1 = (rand(1)*(w_f1_max-w_f1_min))+w_f1_min;
    w_f2 = (rand(1)*(w_f2_max-w_f2_min))+w_f2_min;
    w_phase = (rand(1)*(2*pi-0))+0;
    w_signal(i,:) = w_static_load + w_f1_amp*sin(2*pi*w_f1*t + w_phase) + w_f2_amp*sin(2*pi*w_f2*t + w_phase);
end

%% Simulations

% sensor measure time
s_time = zeros(1,s_qtty);      
for i = 1:s_qtty
    s_time(1,i) = s_pos(i)/v_speed;
end

% sensor measure index
s_idx_ini = zeros(1,s_qtty);
for i = 1:s_qtty
    s_idx_ini(i) = find( abs(t-s_time(i)) < (w_time_res*0.99) );
end

% sensor signal
%s_w_signal = zeros(n,t_size);
s_w_signal = zeros(n,s_qtty); 
s_w_avg = zeros(1,n);

s_linear  = zeros(1,n);
s_pchip   = zeros(1,n);
s_v5cubic = zeros(1,n);
s_makima  = zeros(1,n);
s_spline  = zeros(1,n);

for k = 1:n
    for i = 1:s_qtty
        %s_w_signal(k,s_idx_ini(i)) = w_signal(k,s_idx_ini(i));
        s_w_signal(k,i) = w_signal(k,s_idx_ini(i));
    end
    
    x = interp1(s_time,s_w_signal(k,:),t,'linear');
    x = x(s_idx_ini(1):s_idx_ini(end));
    s_linear(k) = mean(x);
    
    x = interp1(s_time,s_w_signal(k,:),t,'pchip');
    x = x(s_idx_ini(1):s_idx_ini(end));
    s_pchip(k) = mean(x);
    
    x = interp1(s_time,s_w_signal(k,:),t,'v5cubic');
    x = x(s_idx_ini(1):s_idx_ini(end));
    s_v5cubic(k) = mean(x);
    
    x = interp1(s_time,s_w_signal(k,:),t,'makima');
    x = x(s_idx_ini(1):s_idx_ini(end));
    s_makima(k) = mean(x);
    
    x = interp1(s_time,s_w_signal(k,:),t,'spline');
    x = x(s_idx_ini(1):s_idx_ini(end));
    s_spline(k) = mean(x);
    
    %s_linear(k)  = mean(interp1(s_time,s_w_signal(k,:),t,'linear'));
    %s_pchip(k)  = mean(interp1(s_time,s_w_signal(k,:),t,'pchip'));
    %s_v5cubic(k)  = mean(interp1(s_time,s_w_signal(k,:),t,'v5cubic'));
    %s_makima(k)  = mean(interp1(s_time,s_w_signal(k,:),t,'makima'));
    %s_spline(k)  = mean(interp1(s_time,s_w_signal(k,:),t,'spline'));
    
    % Mean value estimator
    s_w_avg(k) = mean(s_w_signal(k,:));
end

%% Outputs
error = zeros(1,n);
error_sig = zeros(1,n);

error_linear  = zeros(1,n);
error_pchip   = zeros(1,n);
error_v5cubic = zeros(1,n);
error_makima  = zeros(1,n);
error_spline  = zeros(1,n);

for i = 1:n
    error(i) = (s_w_avg(i) - w_static_load)*100/w_static_load;
    error_sig(i) = (mean(w_signal(i,:)) - w_static_load)*100/w_static_load;
    
    error_linear(i) = (s_linear(i) - w_static_load)*100/w_static_load;
    error_pchip(i) = (s_pchip(i) - w_static_load)*100/w_static_load;
    error_v5cubic(i) = (s_v5cubic(i) - w_static_load)*100/w_static_load;
    error_makima(i) = (s_makima(i) - w_static_load)*100/w_static_load;
    error_spline(i) = (s_spline(i) - w_static_load)*100/w_static_load;
    
end

x = interp1(s_time,s_w_signal(1,:),t,'pchip');
x = x(s_idx_ini(1):s_idx_ini(end));
plot(t(s_idx_ini(1):s_idx_ini(end)),x)
hold on
x = interp1(s_time,s_w_signal(1,:),t,'makima');
x = x(s_idx_ini(1):s_idx_ini(end));
plot(t(s_idx_ini(1):s_idx_ini(end)),x)
hold on
x = interp1(s_time,s_w_signal(1,:),t,'spline');
x = x(s_idx_ini(1):s_idx_ini(end));
plot(t(s_idx_ini(1):s_idx_ini(end)),x)
hold on
plot(t,w_signal(1,:));

STR = ['Vehicle speed = ',num2str(v_speed*3.6),' km/h, number of runs = ',num2str(n), ', number of sensors = ',num2str(s_qtty)];
disp(STR)
STR = ['MININUM ERROR) mean: ',num2str(min(error),3),' / signal: ',num2str(min(error_sig),3),' / pchip: ',num2str(min(error_pchip),3),' / makima: ',num2str(min(error_makima),3),' / spline: ',num2str(min(error_spline),3)];
disp(STR)
STR = ['MAXIMUM ERROR) mean: ',num2str(max(error),3),' / signal: ',num2str(max(error_sig),3),' / pchip: ',num2str(max(error_pchip),3),' / makima: ',num2str(max(error_makima),3),' / spline: ',num2str(max(error_spline),3)];
disp(STR)
STR = ['MEAN ERROR) mean: ',num2str(mean(error),3),' / signal: ',num2str(mean(error_sig),3),' / pchip: ',num2str(mean(error_pchip),3),' / makima: ',num2str(mean(error_makima),3),' / spline: ',num2str(mean(error_spline),3)];
disp(STR)