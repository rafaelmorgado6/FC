%% aline 2.2
clear all
close all
clc

k = 2;
m =1.5;

a = k/m;
w =sqrt(a);

T = (2*pi())/w;

alfa = -0.2;
t0 = 0;
tf = 10*T;
x0 = 1.9;
v0 = 0;

%tolerancia
reltol = 3e-14;
abstol_1 =1e-13;
abstol_2 =1e-13;

options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t,solucao] = ode45(@f,[t0 tf],[x0 v0],options,k,m,alfa);

%soluçõesk
x = solucao(:,1);
v = solucao(:,2);

plot(t,x);