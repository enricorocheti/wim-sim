clc
close all
clear

t_end = 1;
t = [0:0.00001:t_end-0.00001]; % segundos, resolução de 0.01 ms
f1 = 2; % 1.5 Hz
f2 = 10;  % 10 Hz

static_load = 50;   % kg
x = static_load + sin(2*pi*f1*t) + sin(2*pi*f2*t); % 1 eixo do caminhão

v = 100/3.6;        % km/h / m/s
sensor_length = 7;   % cm
t_sensor = sensor_length/(v*100); % ms

dist = v*t_end; % distância que o veículo percorre durante a duração desse sinal simulado


f_sample = 100000; % 100000 samples por segundo, 1 sample a cada 0.01 ms
n_sample = t_sensor/(1/f_sample); % número de amostras

pos_s1 = 5;  % s1 posicionado em 5 metros
pos_s2 = 10; % s2 posicionado em 10 metros

t_s1_ini = pos_s1/v;
t_s2_ini = pos_s2/v;

idx_s1 = find(abs(t-t_s1_ini)<0.000001);
idx_s2 = find(abs(t-t_s2_ini)<0.000001);

%t_s1 = t(idx_s1:idx_s1+n_sample);
%t_s2 = t(idx_s2:idx_s2+n_sample);

t_size = size(t);
t_size = t_size(2);
x_s1 = zeros(1,t_size);
x_s2 = zeros(1,t_size);

for i=idx_s1:idx_s1+n_sample
    x_s1(i) = x(i);
end
for i=idx_s2:idx_s2+n_sample
    x_s2(i) = x(i);
end

figure(1);
%plot(t,x);
plot(t(15000:40000),x(15000:40000),"b");
hold on
stem(t(15000:40000),x_s1(15000:40000),"g","filled");
hold on
stem(t(15000:40000),x_s2(15000:40000),"y","filled");
ylim([0,static_load+5]);
xlabel('t [s]');
ylabel('Carga dinâmica [kg]');
%yline(static_load,'-','Carga estática');
line ([0.15 0.4], [static_load static_load], "linestyle", "--", "color", "r");


x_sensor = x_s1 + x_s2;
X = fft(x);
figure(2);
plot(abs(X));
