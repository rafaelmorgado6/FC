clear all, clc

h=0.01;
x0=1;
v0=0;
K=16;
m=1;
tf=10;

reltol = 3e-14;
abstol_1 = 1e-13;
abstol_2 = 1e-13;

options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t,solucao]=ode45(@f,0:h:tf,[x0,v0],options,m,K);

tt=t';
xx=solucao(:,1)';
vv=solucao(:,2)';
figure(1);
plot(tt,xx)
figure(2);
plot(xx,vv)

Ecinetica=1/2.*m*vv.^2;
Epotencial=1/2.*K.*xx.^2;
Emecanica=Ecinetica+Epotencial;

figure(3)
plot(tt,Emecanica);
title('Emec em função do número de iterações');