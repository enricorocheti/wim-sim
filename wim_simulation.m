clc
close all
clear
rng(10)                                         % fix seed for rand function

%% Configurations
config = jsondecode(fileread('config.json'));	% config file

v_speeds = config.vehicle_speed./3.6;           % vehicle speed (km/h -> m/s) 
v_qtty   = length(config.vehicles);
s_qtty   = config.sensor_qtty;                  % number of sensors
s_dist   = config.sensor_distance;              % distance between sensors (meters)
s_rsd    = config.sensor_rsd;                   % sensors' relative standard deviation
n_sim    = config.number_of_runs;               % number of simulations 
v_max_ax = 7;                                 	% hard-coded value of the maximum number of axles

for v_speed = v_speeds(:).'    
    
%% System parameters and layout
s_pos = zeros(1,s_qtty);        % sensor position (meters)
for i = 1:s_qtty
    s_pos(1,i) = i * s_dist;    % TODO: s_pos starting on zero? now it's not
end

% f1 = (5 - 1)/2 = 3 Hz
% f2 = (15 - 8)/2 = 11.5 Hz
s_design_1 = 2*(s_qtty - 1)*v_speed/(3*s_qtty^2);
s_design_2 = (v_speed/(2*s_qtty))*(1/3 + (s_qtty-1)/11.5);

%% Vehicle load signal parameters
w_time_end = s_pos(end)/v_speed;                % time that the vehicle travel through all sensors
w_time_end = w_time_end + 0.1;                  % adding some time as a margin

w_time_res = 0.001;                             % time resolution (seconds)
t = (0:w_time_res:w_time_end-w_time_res);       % time vector
t_size = size(t,2);                             % time vector size

% dynamic loads frequencies range (Hz)
f1_min = 1;               
f1_max = 5;
f2_min = 8;
f2_max = 15;

% dynamic load amplitudes (tons) based on vehicle speed linearization
w1 = 0.004167*(v_speed*3.6) + 0.015000; 
w2 = 0.001042*(v_speed*3.6) + 0.003750;

% vector containing signals from each vehicle's axles
w_signal = zeros(n_sim,v_qtty,v_max_ax,t_size);

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
            %plot(t,w_signal_axle);
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
s_w_signal = zeros(n_sim, v_qtty, v_max_ax, s_qtty);

% estimators
axle_mean       = zeros(n_sim, v_qtty, v_max_ax);
axle_MLE        = zeros(n_sim, v_qtty, v_max_ax);
axle_sr_pchip   = zeros(n_sim, v_qtty, v_max_ax);
axle_sr_makima  = zeros(n_sim, v_qtty, v_max_ax);
axle_sr_spline  = zeros(n_sim, v_qtty, v_max_ax);

gvw_mean        = zeros(n_sim, v_qtty);
gvw_MLE         = zeros(n_sim, v_qtty); 
gvw_sr_pchip    = zeros(n_sim, v_qtty);
gvw_sr_makima   = zeros(n_sim, v_qtty);
gvw_sr_spline   = zeros(n_sim, v_qtty);

for i = 1:n_sim
    for j = 1:v_qtty
        for k = 1:config.vehicles(j).axle_qtty
            numer_MLE = 0;
            denom_MLE = 0;
            
            for l = 1:s_qtty
                sample = w_signal(i, j, k, s_idx_ini(l));
                current_std = sample * s_rsd(l);
                current_var = current_std^2;
                
                % apply normally distributed noise using the current standard deviation
                s_w_signal(i,j,k,l) = normrnd(sample, current_std);
                
                % MLE computation
                numer_MLE = numer_MLE + (s_w_signal(i,j,k,l)/current_var);
                denom_MLE = denom_MLE + (1/current_var);
                
                % plot sensor samples
                %stem(s_time(l),s_w_signal(i,j,k,l));
                %hold on
            end
            
            % mean value
            axle_mean(i,j,k) = mean(s_w_signal(i,j,k,:));
            gvw_mean(i,j) = gvw_mean(i,j) + axle_mean(i,j,k);
            
            % MLE
            axle_MLE(i,j,k) = numer_MLE / denom_MLE;
            gvw_MLE(i,j) = gvw_MLE(i,j) + axle_MLE(i,j,k);
                
            % signal reconstruction
            signal_2d = reshape(s_w_signal(i, j, k, :), 1, []);
            
            x = interp1(s_time,signal_2d,t,'pchip');
            axle_sr_pchip(i,j,k) = mean(x(s_idx_ini(1):s_idx_ini(end))); 
            gvw_sr_pchip(i,j) = gvw_sr_pchip(i,j) + axle_sr_pchip(i,j,k);
            
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

% static gvw
v_static_gvw = zeros(v_qtty);
for i = 1:v_qtty
    for j = 1:config.vehicles(i).axle_qtty
        v_static_gvw(i) = v_static_gvw(i) + config.vehicles(i).axle_st_load(j);
    end
end

err_axl_mv      = zeros(n_sim,v_qtty,v_max_ax);
err_axl_MLE   	= zeros(n_sim,v_qtty,v_max_ax);
err_axl_pchip   = zeros(n_sim,v_qtty,v_max_ax);
err_axl_makima  = zeros(n_sim,v_qtty,v_max_ax);
err_axl_spline  = zeros(n_sim,v_qtty,v_max_ax);

err_gvw_mv      = zeros(n_sim,v_qtty);
err_gvw_MLE   	= zeros(n_sim,v_qtty);
err_gvw_pchip   = zeros(n_sim,v_qtty);
err_gvw_makima  = zeros(n_sim,v_qtty);
err_gvw_spline  = zeros(n_sim,v_qtty);

for i = 1:n_sim
    for j = 1:v_qtty
        for k = 1:config.vehicles(j).axle_qtty
            % axle error
            err_axl_mv(i,j,k)      = (axle_mean(i,j,k)       - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_MLE(i,j,k)     = (axle_MLE(i,j,k)        - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            
            err_axl_pchip(i,j,k)   = (axle_sr_pchip(i,j,k)   - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_makima(i,j,k)  = (axle_sr_makima(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
            err_axl_spline(i,j,k)  = (axle_sr_spline(i,j,k)  - config.vehicles(j).axle_st_load(k)) * 100 / config.vehicles(j).axle_st_load(k);
        end
        
        % gvw error
        err_gvw_mv(i,j)      = (gvw_mean(i,j)       - v_static_gvw(j)) * 100 / v_static_gvw(j);
        err_gvw_MLE(i,j)     = (gvw_MLE(i,j)        - v_static_gvw(j)) * 100 / v_static_gvw(j);
        err_gvw_pchip(i,j)   = (gvw_sr_pchip(i,j)   - v_static_gvw(j)) * 100 / v_static_gvw(j);
        err_gvw_makima(i,j)  = (gvw_sr_makima(i,j)  - v_static_gvw(j)) * 100 / v_static_gvw(j);
        err_gvw_spline(i,j)  = (gvw_sr_spline(i,j)  - v_static_gvw(j)) * 100 / v_static_gvw(j);
    end
end

%% CSV output
% axle estimation statistics
csvName = ['outputs/axle_output_s',num2str(s_qtty),'_d',num2str(s_dist,2),'.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');   
    fprintf(csvFile, 'speed,');
    fprintf(csvFile, 'mean_mv,mean_MLE,mean_pchip,mean_makima,mean_spline,');
    fprintf(csvFile, 'max_mv,max_MLE,max_pchip,max_makima,max_spline,');
    fprintf(csvFile, 'std_mv,std_MLE,std_pchip,std_makima,std_spline\n');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '%.0f,', v_speed * 3.6);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,', mean(err_axl_mv(:)), mean(err_axl_MLE(:)), mean(err_axl_pchip(:)), mean(err_axl_makima(:)), mean(err_axl_spline(:)));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,', max(abs(err_axl_mv(:))), max(abs(err_axl_MLE(:))), max(abs(err_axl_pchip(:))), max(abs(err_axl_makima(:))), max(abs(err_axl_spline(:))));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f\n', std(err_axl_mv(:)), std(err_axl_MLE(:)), std(err_axl_pchip(:)), std(err_axl_makima(:)), std(err_axl_spline(:)));
fclose(csvFile);

% GVW estimation statistics
csvName = ['outputs/gvw_output_s',num2str(s_qtty),'_d',num2str(s_dist,2),'.csv'];
if ~exist(csvName,'file')
    csvFile = fopen(csvName, 'w');   
    fprintf(csvFile, 'speed,');
    fprintf(csvFile, 'mean_mv,mean_MLE,mean_pchip,mean_makima,mean_spline,');
    fprintf(csvFile, 'max_mv,max_MLE,max_pchip,max_makima,max_spline,');
    fprintf(csvFile, 'std_mv,std_MLE,std_pchip,std_makima,std_spline\n');
else
    csvFile = fopen(csvName, 'a');
end

fprintf(csvFile, '%.0f,', v_speed * 3.6);
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,', mean(err_gvw_mv(:)), mean(err_gvw_MLE(:)), mean(err_gvw_pchip(:)), mean(err_gvw_makima(:)), mean(err_gvw_spline(:)));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f,', max(abs(err_gvw_mv(:))), max(abs(err_gvw_MLE(:))), max(abs(err_gvw_pchip(:))), max(abs(err_gvw_makima(:))), max(abs(err_gvw_spline(:))));
fprintf(csvFile, '%.3f,%.3f,%.3f,%.3f,%.3f\n', std(err_gvw_mv(:)), std(err_gvw_MLE(:)), std(err_gvw_pchip(:)), std(err_gvw_makima(:)), std(err_gvw_spline(:)));
fclose(csvFile);

end