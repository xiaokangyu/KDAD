% function ZHAI(tm,dt)
%tm为仿真时间，dt为仿真时间步长
clc;
clear;
close all;%关闭所有的Figure窗口

tm=1.0;
dt=0.00001;
tspan = [0 1];
y0 = zeros(32,1);
[t,dx] = ode45('MyBearingFun',tspan,y0);
plot(dx(:,17))
for i = 17:32
figure
f = abs(fft(dx(:,i)))/length(dx(:,6));
f(1)=0;
plot(f(1:2500))
end
