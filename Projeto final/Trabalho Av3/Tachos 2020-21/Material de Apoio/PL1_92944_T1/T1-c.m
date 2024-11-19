%Ema Fadiga - 92944 - PL1

%% c - Amplitude e Período em função de mu
clear all 
close all
clc

m=1;
K=1;
alpha=0.2;
y(1)=1.5;
v(1)=0;
tf=150;
F0=0;
w0=1;
% w0=2;  % Optar por um dos w0

dt=0.1;
t=0:dt:tf;
Nt=numel(t);

muu=0:0.05:0.8;
Nu=numel(muu);

A=zeros(Nu,1);
P=zeros(Nu,1);


for k=1:Nu
    mu=muu(k);
    options=odeset('RelTol', 1e-13, 'AbsTol', [1e-13 1e-13]);
    [t, solution]=ode45(@ODE, t, [y(1) v(1)], options,m,K,alpha,mu,w0,F0);
  
    y=solution(:,1);
    v=solution(:,2);
    Periodos=0;
    imax=0;
    for i =2:Nt-1
        if and (y(i-1)<y(i),y(i)>y(i+1))
            imax=imax+1;
            aux=lagr(t(i-1:i+1),y(i-1:i+1));
            tmax(imax)=aux(1); 
            ymax(imax)=aux(2);
        end
    end
    ym=mean(ymax);
    A(k)= ym;
    for j=2:numel(tmax)
       Periodos(j-1)= tmax(j)-tmax(j-1);
    end
    P(k)= mean(Periodos);
    clear ymax
    clear tmax
end

figure(1)
plot(muu,A,'b*')
title('Variação da Amplitude')
xlabel('mu')
ylabel('Amplitude')

hold on

figure(2)
plot(muu,P,'r*')
title('Variação do Período')
xlabel('mu')
ylabel('Período')
