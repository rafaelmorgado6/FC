clc;
clear all;
close all; 

t0 = 0;
tf = 1.6;
h = 0.2;
t = t0:h:tf;
N = length(t);
v0 = 0;
m = 0.150;
g = 9.8;
v = zeros(N,1);
v(1) = v0;

% Método de Euler

k = 1;
while k < N
 
    v(k+1) = v(k) + (-g)*h ;
    k=k+1;
end

% solução analitica

vz = v0 - g*t;

plot(t,v,t,vz,'*')
grid on
title('Comparação do método de Euler com a solução analitica')
xlabel('tempo (s)')
ylabel('velocidade (m/s)')
legend('Método de Euler','Solução Analítica')