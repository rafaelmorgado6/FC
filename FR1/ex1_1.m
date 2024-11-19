%%                                  Ficha de Revisões 1
clear all;close all; clc;
% Exercício 1 -> Crank-Nicolson
% variáveis
m=1;
K=1;
alfa=-0.1;

%condições iniciais
x(1)=1;
vx(1)=1;

%Método de Crank-Nicolson
h=0.01;
T=2*pi*sqrt(m/K);
tf=10;
t=0:h:tf;
N=numel(t);

%opções pedidas no guião
const = [h/2, K*h/(2*m), 2*alfa];
options = optimset('Display','off','Tolx',1e-10,'TolFun',1e-10);

imax=0;

for k=1:N-1
    %chamar função fcr
    func = @(xv) fcr(xv,x(k),vx(k),const);
    % aplicar método
    xv0 = [x(k),vx(k)];
    aux = fsolve(func,xv0,options);
    x(k+1) = aux(1);
    vx(k+1) = aux(2);
end   

% calcular amplitude e periodo
for i=2:N-1
    if and (x(i+1)-x(i)<=0, x(i)-x(i-1)>=0)
        imax = imax + 1;
        aprox = lagr(t(i-1:i+1),x(i-1:i+1));
        tmax(imax)=aprox(1);
        xmax(imax)=aprox(2);
    end
end

p=polyfit(1:imax,tmax,1)
Periodo=p(1)
Amplitude=mean(xmax)
    
figure(1)
plot(t,x)

figure(2)
plot(t,vx)

figure(3)
plot(x,vx)