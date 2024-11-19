%% Folha Revisões 2
clear all;close all; clc;
v = 1.7187;
B = 18;
h = 0.01;
x = 0;
int = 0:h:7;
T = zeros(1,length(int));
dT = zeros(1,length(int));

T(1) = 10^-4;
dT(1) = 10^-4*v;

for i=1:length(int)
    dT(i+1)= dT(i) + ((v-B*exp(-1/T(i)))*dT(i)-B*exp(-1/T(i))*(1-v*T(i)))*h;
    T(i+1) = T(i) + dT(i+1)*h; 
end

figure()
plot(int,T(1:end-1))
figure()
plot(int,dT(1:end-1))