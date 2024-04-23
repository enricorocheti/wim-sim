clc
close all
clear
rng('shuffle')                                  % generate new seed for rand function

%% Configurations
config = jsondecode(fileread('config.json'));	% config file

v_speeds = config.vehicle_speed./3.6;           % vehicle speed (km/h -> m/s) 
v_qtty   = length(config.vehicles);
s_qtty   = config.sensor_qtty;                  % number of sensors
s_dist   = config.sensor_distance;              % distance between sensors (meters)
s_rsd    = config.sensor_rsd;                   % sensors' relative standard deviation
n_sim    = config.number_of_runs;               % number of simulations 

for v_speed = v_speeds(:).'    
    
%% System parameters and layout
s_pos = zeros(1,s_qtty);        % sensor position (meters)
for i = 1:s_qtty
    s_pos(1,i) = i * s_dist;    % TODO: s_pos starting on zero? now it aint
end

% f1 = (5 - 1)/2 = 3 Hz
% f2 = (15 - 8)/2 = 11.5 Hz
s_design_1 = 2*(s_qtty - 1)*v_speed/(3*s_qtty^2);
s_design_2 = (v_speed/(2*s_qtty))*(1/3 + (s_qtty-1)/11.5);

%% Vehicle load signal parameters
w_time_end = s_pos(end)/v_speed;    % time that the vehicle travel through all sensors
w_time_end = w_time_end + 0.1;      % adding some time as a margin

w_time_res = 0.001;                             % time resolution (seconds)
t = (0:w_time_res:w_time_end-w_time_res);       % time vector
t_size = size(t,2);                             % time vector size

% dynamic load frequencies range (Hz)
f1_min = 1;               
f1_max = 5;
f2_min = 8;
f2_max = 15;

% dynamic load amplitudes (tons) based on vehicle speed linearization
w1 = 0.004167*(v_speed*3.6) + 0.015000; 
w2 = 0.001042*(v_speed*3.6) + 0.003750;

% vector containing signals from each vehicle's axles
w_signal = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty),t_size);

for i = 1:n_sim
    for j = 1:v_qtty
        f1 = (rand(1)*(f1_max-f1_min)) + f1_min;
        f2 = (rand(1)*(f2_max-f2_min)) + f2_min;
        phase = (rand(1)*(2*pi-0)) + 0;

        axle_qtty = config.vehicles(j).axle_qtty;
        for k = 1:axle_qtty
            st_load = config.vehicles(j).axle_st_load(k);
            if k == 1
                w_signal_axle = st_load*(1 + w1*sin(2*pi*f1*t + phase) + w2*sin(2*pi*f2*t + phase));
            else
                deltaD = abs(config.vehicles(j).axle_dist(k) - config.vehicles(j).axle_dist(k-1));
                deltaT = deltaD / v_speed;
                w_signal_axle = st_load*(1 + w1*sin(2*pi*f1*(t-deltaT) + phase) + w2*sin(2*pi*f2*(t-deltaT) + phase));
            end
            w_signal(i, j, k, :) = w_signal_axle;
            
            % plot axle signal
            plot(t,w_signal_axle);
            hold on
        end
    end
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
    x = find( abs(t-s_time(i)) < (w_time_res*0.99) );
    s_idx_ini(i) = x(1);
end

% sensor signal
s_w_signal = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty), s_qtty);

% estimators
axle_mean   = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
axle_sr_linear  = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
axle_sr_pchip   = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
axle_sr_v5cubic = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
axle_sr_makima  = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
axle_sr_spline  = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));

gvw_mean = zeros(n_sim, v_qtty);
gvw_sr_linear = zeros(n_sim, v_qtty);
gvw_sr_pchip = zeros(n_sim, v_qtty);
gvw_sr_v5cubic = zeros(n_sim, v_qtty);
gvw_sr_makima = zeros(n_sim, v_qtty);
gvw_sr_spline = zeros(n_sim, v_qtty);

for i = 1:n_sim
    for j = 1:v_qtty
        for k = 1:config.vehicles(j).axle_qtty
            for l = 1:s_qtty
                sample = w_signal(i, j, k, s_idx_ini(l));
                current_std = sample * s_rsd(l);
                current_var = current_std^2;
                
                % apply normally distributed noise using the current standard deviation
                s_w_signal(i,j,k,l) = normrnd(sample, current_std);
                
                % plot sensor samples
                stem(s_time(l),s_w_signal(i,j,k,l));
                hold on
            end
            
            % Mean value
            axle_mean(i,j,k) = mean(s_w_signal(i,j,k,:));
            gvw_mean(i,j) = gvw_mean(i,j) + axle_mean(i,j,k);
                
            % Signal reconstruction
            signal_2d = reshape(s_w_signal(i, j, k, :), 1, []);
            
            x = interp1(s_time,signal_2d,t,'linear');
            axle_sr_linear(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));
            gvw_sr_linear(i,j) = gvw_sr_linear(i,j) + axle_sr_linear(i,j,k);
            
            x = interp1(s_time,signal_2d,t,'pchip');
            axle_sr_pchip(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end))); 
            gvw_sr_pchip(i,j) = gvw_sr_pchip(i,j) + axle_sr_pchip(i,j,k);
            
            x = interp1(s_time,signal_2d,t,'v5cubic');
            axle_sr_v5cubic(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));  
            gvw_sr_v5cubic(i,j) = gvw_sr_v5cubic(i,j) + axle_sr_v5cubic(i,j,k);
            
            x = interp1(s_time,signal_2d,t,'makima');
            axle_sr_makima(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));  
            gvw_sr_makima(i,j) = gvw_sr_makima(i,j) + axle_sr_makima(i,j,k);
            
            x = interp1(s_time,signal_2d,t,'spline');
            axle_sr_spline(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));
            gvw_sr_spline(i,j) = gvw_sr_spline(i,j) + axle_sr_spline(i,j,k);
        end
    end
end

%% Outputs

v_static_gvw = zeros(v_qtty);
for i = 1:v_qtty
    for j = 1:config.vehicles(i).axle_qtty
        v_static_gvw(i) = v_static_gvw + config.vehicles(i).axle_st_load(j);
    end
end

%err_gvw_mv = zeros(n_sim,v_qtty);

err_axl_mv      = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty));
err_axl_linear  = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty));
err_axl_pchip   = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty));
err_axl_v5cubic = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty));
err_axl_makima  = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty));
err_axl_spline  = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty));
for i = 1:n_sim
    for j = 1:v_qtty
        for k = 1:config.vehicles(j).axle_qtty
            err_axl_mv(i,j,k)      = (axle_mean(i,j,k)   - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            
            err_axl_linear(i,j,k)  = (axle_sr_linear(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_pchip(i,j,k)   = (axle_sr_pchip(i,j,k)   - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_v5cubic(i,j,k) = (axle_sr_v5cubic(i,j,k) - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_makima(i,j,k)  = (axle_sr_makima(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_spline(i,j,k)  = (axle_sr_spline(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
        end
        %err_gvw_mv(i,j) = (v_gvw(i,j) - v_static_gvw(j)) * 100 / v_static_gvw(j);
    end
end

%% Terminal output
STR = ['Speed = ',num2str(v_speed),' m/s, Static load'];
disp(STR)
STR = ['f1 = ',num2str(f1),' Hz, A1 = ',num2str(st_load * w1),' Kg'];
disp(STR)
STR = ['f2 = ',num2str(f2),' Hz, A2 = ',num2str(st_load * w2),' Kg'];
disp(STR)

%% CSV output

csvName = ['output_s',num2str(s_qtty),'_n',num2str(n_sim),'.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');
    fprintf(csvFile, 'speed,mean_mv,mean_pchip,mean_makima,mean_spline,std_mv,std_pchip,std_makima,std_spline,min_mv,min_pchip,min_makima,min_spline,max_mv,max_pchip,max_makima,max_spline\n');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '%.0f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n', v_speed*3.6,mean(err_axl_mv),mean(err_axl_pchip),mean(err_axl_makima),mean(err_axl_spline),std(err_axl_mv),std(err_axl_pchip),std(err_axl_makima),std(err_axl_spline),min(err_axl_mv),min(err_axl_pchip),min(err_axl_makima),min(err_axl_spline),max(err_axl_mv),max(err_axl_pchip), max(err_axl_makima),max(err_axl_spline));
fclose(csvFile);

end

return

fprintf(csvFile, 'sim_index,vehicle_index,axle_index,mean_value,err_mean_value\n');
for i = 1:n_sim
    for j = 1:v_qtty
        for k = 1:config.vehicles(j).axle_qtty
            fprintf(csvFile, '%d,%d,%d,%.3f,%.3f\n', i, j, k, axle_mean(i,j,k), err_axl_mv(i,j,k));
        end
    end
end
fclose(csvFile);

return