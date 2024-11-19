%Ema Fadiga - 92944 - PL1

%% b- ODE45 
clear all 
close all

m=1;
K=1;
alpha=0.2;
mu=0.8;
y(1)=1.5;
v(1)=0;
tf=150;
% w0=1;
w0=2;  % Optar por um dos w0
% F0=0; 
F0=0.18;  % Optar por um dos F0

options=odeset('RelTol', 1e-13, 'AbsTol', [1e-13 1e-13]);
[t, solution]=ode45(@ODE,[0 tf], [y(1) v(1)], options,m,K,alpha,mu,w0,F0);
 
y=solution(:,1);
v=solution(:,2);% mostra a solução com o caos

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
ylabel('v')% mostra a solução com o caos

hold on
figure(3)

N=numel(t)%t é dado pela ode45
N1=round(0.95*N);%é 95% do ciclo inicial 
figure(3)
plot(y(N1:N),v(N1:N),'k')%representa o ciclo limite
    
title('Espaço de fases')
xlabel('y')
ylabel('v')% mostra a solução com o caos

