clc
close all
clear

config = jsondecode(fileread('config.json'));   % config file

v_speeds = config.vehicle_speed./3.6;           % vehicle speed (km/h -> m/s) 
v_qtty = length(config.vehicles);
s_qtty = config.sensor_qtty;                    % number of sensors
s_dist = config.sensor_distance;                % distance between sensors (meters)
n_sim = config.number_of_runs;                  % number of simulations 

for v_speed = v_speeds(:).'

%% System parameters and layout
s_pos = zeros(1,s_qtty);        % sensor position (meters)
for i = 1:s_qtty
    s_pos(1,i) = i * s_dist;    % TODO: s_pos starting on zero? now it aint
end

%% Weight signal parameters
w_time_end = s_pos(end)/v_speed;    % time to the vehicle travel through all sensors
w_time_end = w_time_end + 0.1;      % adding time as a margin

w_time_res = 0.001;                           % time resolution (seconds)
t = (0:w_time_res:w_time_end-w_time_res);       % time vector
t_size = size(t,2);                             % time vector size

% f1 and f2 range (Hz)
f1_min = 1;               
f1_max = 5;
f2_min = 8;
f2_max = 15;

% TODO: amplitude should to depend on vehicle speed
w_f1_amp = 1;     	% first dynamic load amplitude (Kg)    
w_f2_amp = 0.5;     % second dynamic load amplitude (Kg) 

% vector containing signals from each vehicle's axles
w_signal = zeros(n_sim,v_qtty,max(config.vehicles.axle_qtty),t_size);

for i = 1:n_sim
    for j = 1:v_qtty
        w_f1 = (rand(1)*(f1_max-f1_min)) + f1_min;    % TODO: seed deve ser sempre aleatório
        w_f2 = (rand(1)*(f2_max-f2_min)) + f2_min;    % TODO: seed deve ser sempre aleatório
        w_phase = (rand(1)*(2*pi-0)) + 0;                   % TODO: seed deve ser sempre aleatório

        axle_qtty = config.vehicles(j).axle_qtty;
        
        for k = 1:axle_qtty
            w_signal_axle = config.vehicles(j).axle_st_load(k) + w_f1_amp*sin(2*pi*w_f1*t + w_phase) + w_f2_amp*sin(2*pi*w_f2*t + w_phase);
            w_signal(i, j, k, :) = w_signal_axle;
            %plot(w_signal_axle)
            %hold on
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

% mean value estimator
s_w_avg   = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));

% signal reconstruction estimators
s_linear  = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
s_pchip   = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
s_v5cubic = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
s_makima  = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));
s_spline  = zeros(n_sim, v_qtty, max(config.vehicles.axle_qtty));

v_gvw = zeros(n_sim, v_qtty);

for i = 1:n_sim
    for j = 1:v_qtty
        for k = 1:config.vehicles(j).axle_qtty
            for l = 1:s_qtty
                % TODO: add time offset between axles
                s_w_signal(i,j,k,l) = w_signal(i, j, k, s_idx_ini(l));
            end
            
            % Gross vehicle weight
            %v_gvw(i,j) = v_gvw(i,j) + config.vehicles(j).axle_st_load(k);
            
            % Mean value
            s_w_avg(i,j,k) = mean(s_w_signal(i,j,k,:));
                
            % Signal reconstruction
            signal_2d = reshape(s_w_signal(i, j, k, :), 1, []);
            
            x = interp1(s_time,signal_2d,t,'linear');
            s_linear(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));
            x = interp1(s_time,signal_2d,t,'pchip');
            s_pchip(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end))); 
            x = interp1(s_time,signal_2d,t,'v5cubic');
            s_v5cubic(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));  
            x = interp1(s_time,signal_2d,t,'makima');
            s_makima(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));  
            x = interp1(s_time,signal_2d,t,'spline');
            s_spline(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end)));
        end
    end
end

%% Outputs

%v_static_gvw = zeros(v_qtty);
%for i = 1:v_qtty
%    for j = 1:config.vehicles(i).axle_qtty
%        v_static_gvw(i) = v_static_gvw + config.vehicles(i).axle_st_load(j);
%    end
%end

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
            err_axl_mv(i,j,k)      = (s_w_avg(i,j,k)   - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            
            err_axl_linear(i,j,k)  = (s_linear(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_pchip(i,j,k)   = (s_pchip(i,j,k)   - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_v5cubic(i,j,k) = (s_v5cubic(i,j,k) - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_makima(i,j,k)  = (s_makima(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_spline(i,j,k)  = (s_spline(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
        end
        %err_gvw_mv(i,j) = (v_gvw(i,j) - v_static_gvw(j)) * 100 / v_static_gvw(j);
    end
end

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
            fprintf(csvFile, '%d,%d,%d,%.3f,%.3f\n', i, j, k, s_w_avg(i,j,k), err_axl_mv(i,j,k));
        end
    end
end
fclose(csvFile);

return