%Ema Fadiga - 92944 - PL1

%% d - RK4
clear all 
close all

m=1;
K=1;
alpha=0.2;
mu=0.8;
y(1)=1.5;
v(1)=0;
tf=150;
w0=1;
% w0=2;  % Optar por um dos w0
F0=0; 

dt=0.001;
tf=150;
t=0:dt:tf;
N=numel(t);

fy=@(t,y,v) v;
fv=@(t,y,v) (-K*y*(1+alpha*y^2) +mu*cos(v)*v + F0*cos(w0*t))/m;

for i=1:N-1
    r1v=fv(t(i),y(i),v(i));
    r1y=fy(t(i),y(i),v(i));
    
    r2v=fv( t(i)+dt/2, y(i)+r1y*dt/2 , v(i)+r1v*dt/2);
    r2y=fy( t(i)+dt/2, y(i)+r1y*dt/2 , v(i)+r1v*dt/2);
    
    r3v=fv( t(i)+dt/2, y(i)+r2y*dt/2, v(i)+r2v*dt/2);
    r3y=fy( t(i)+dt/2, y(i)+r2y*dt/2, v(i)+r2v*dt/2);
    
    r4v=fv( t(i)+dt, y(i)+r3y*dt, v(i)+r3v*dt);
    r4y=fy( t(i)+dt, y(i)+r3y*dt, v(i)+r3v*dt);
    
    y(i+1)=y(i)+dt*(r1y+2*r2y+2*r3y+r4y)/6;
    v(i+1)=v(i)+dt*(r1v+2*r2v+2*r3v+r4v)/6;
end


figure(1)
plot(t,y)
title('y(t)')
xlabel('t')
ylabel('y')

hold on

figure(2)
plot(y,v,'r')
title('Espaço de fases')
xlabel('y')
ylabel('v')

