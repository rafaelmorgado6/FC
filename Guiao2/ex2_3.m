%%                                      Trabalho Prático 2
close all; clear all; clc;
%% Exercício 2
% alinea a
close all; clear all; clc;
x(1)=0.47;
v(1)=8.2;
 m=1;
% Método de Euler-Cromer
h= 0.0001;
tf = 50;
t=0:h:tf;
N=numel(t);

for i=1:N-1
    
    f(i)=-k*x(i)/m;
    v(i+1) = v(i)+f(i);
    x(i+1) = x(i)+v(i+1)*h;
end

m=1;
K=1;
% V=1/2*K*x^2*(1+alpha*x^2);
% F=-K*x*(1+2*alpha*x^2);

%% a
alpha=-0.1;
x(1)=1;
v(1)=1;
h=0.0001;

for k=1:1000000
    a(k+1)=-K/m*(x(k)+2*alpha*x(k)^3);
    v(k+1)=v(k)-((x(k)^2*(x(k)^2 - 10))/20)*h;
    x(k+1)=x(k)+((x(k)^3*(3*x(k)^2 - 50))/300)*h;
end

t=0:h:k*h;
figure(1)
plot(t,x)
title('X em função de t');

figure(2)
plot(t,v)
title('V em função de t');

figure(3)
plot(x,v)
title('V em função de X');