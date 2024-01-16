clear all
close all
clc

fa  = 100;          % freq. de amostragem
dt  = 1/fa;          % discretização no tempo
t   = [0:dt:1-dt];   % vetor tempo

M = csvread('teste.csv')
M(:,1) = []

figure(1);
t = M(:,1); % tempo
p = M(:,2); % pressão da roda sobre o sensor
t = t(15000:32000);
p = p(15000:32000);
plot(t,p);

N = size(p);
N(:,2) = [];
%P = fft(p.*hanning(N)')/N;

Ts = 0.05;
Fs = 1/Ts;
[y, Ty] = resample(p,t,Fs);
lowpass(y,0.5,Fs)

return

dt = (t(end)-t(1))/N;
df = 1/N/dt;
f = 0:df:(N-1)*df
%Xp = P(1:(size(p)/2));
%Xp(2:end-1) = 2*Xp(2:end-1);
%figure(2);
%plot(abs(Xp));

figure(2);
%plot(f(1:N/2),abs(P(1:N/2)));
title('Espectro original do sinal');
xlabel('f [Hz]');
    
