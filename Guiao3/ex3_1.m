clear all; 
close all; 
clc

m = 1; 
K = 16;
tf = 10; 
A = 1;
W = sqrt(K/m);
h = 0.1;
t = 0:h:tf;

vx = zeros(1,length(t));
x = zeros(1,length(t));
Em = zeros(1,length(t));
x(1) = 1;


% Funções anónimas
fx = @(v) v;
fv = @(x) -K*x/m; 


%% Runge-Kutta 2nd order
for i = 1:(length(t)-1)
    r1v = fv(x(i));
    r1x = fx(vx(i));
    
    r2v = fv(x(i)+ r1x* (h/2));
    r2x = fx(vx(i) + r1v * (h/2));
    
    vx(i+1) = vx(i) + r2v*h;
    x(i+1) = x(i) + r2x*h;
    
    Em(i) = (1/2)*m*abs(vx(i))^2 + (1/2)*K*x(i)^2;
end

%% Solução Analítica
x_a = A*cos(W*t); 
vx_a = A*W*sin(W*t);
Em_a = (1/2)*K*A^2;

figure(1)
subplot(4,1,1)
plot(t,x,'r',t,x_a,'b')
subplot(4,1,2)
plot(t,vx,'y',t,vx_a,'b')
subplot(4,1,3)
plot(x,vx,'b',x,vx_a,'r')
subplot(4,1,4)
plot(t,Em,'g',t,Em_a,'b')

%% Euler
for i = 1:(length(t)-1)
    ax = -(K/m)*x(i);
    vx(i+1) = vx(i) + ax*h;
    x(i+1) = x(i) + vx(i)*h;
    
    Em(i) = (1/2)*m*abs(vx(i))^2 + (1/2)*K*x(i)^2;
end

figure(2)
subplot(4,1,1)
plot(t,x,'r')
subplot(4,1,2)
plot(t,vx,'y')
subplot(4,1,3)
plot(x,vx,'b')
subplot(4,1,4)
plot(t,Em,'g')



