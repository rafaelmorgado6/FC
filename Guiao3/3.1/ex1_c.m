%%  Ex1.c

clear,clc,close all

%% CONSTANTES

h=0.01;
tf=10;
t=0:h:tf;
N=length(t); 

k=16;
m=1;

x=[];
v=[];

x(1)=1;
y(1)=10;

yexact=3*exp(-2.*t); %solucao exata

%% MÉTODO DE RUNGE-KUTTA
    for k=1:N-1
k1=-2*y(k);
 y1=y(k)+k1*h/2;
k2=-2*y1;
 y2=y(k)+k2*h/2;
k3=-2*y2;
 y3=y(k)+k3*h;
k4=-2*y3;
 y(k+1)=y(k)+h/6*(k1+2*k2+2*k3+k4);
end
plot(t,yexact,'*',t,y)