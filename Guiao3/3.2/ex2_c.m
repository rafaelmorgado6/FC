%%  Ex2.c

clear,clc,close all

%% CONSTANTES

h=0.5;  %TROCAR AQUI O h
tf=10;      
t=0:h:tf;

k=16;
m=1;

w=sqrt(k/m);

x=[];
v=[];

x(1)=1;
v(1)=0;

%% MÉTODO DE RUNGE-KUTTA de 4ª ordem
%PARA SIMPLIFICAR SÓ É NECESSÁRIO COLOCAR AS VARIAVEIS DAS EXPRESSOES, NÃO
%É NECESSÁRIO COLOCAR (t,X,V)

fx = @(V) V;        %derivada em ordem ao tempo de x
fv = @(X) -k*X/m;   %derivada em ordem ao tempo de v

for i=1:length(t)-1 
   r1x=fx(t+0*h,x(i),v(i));    %1ª derivada *0
   r1v=fv(x(i));

   r2x=fx(v(i)+h/2*r1v);    %2ª derivada *1/2
   r2v=fv(x(i)+h/2*r1x);
    
   r3x=fx(v(i)+h/2*r2v);    %3ª derivada *1/2
   r3v=fv(x(i)+h/2*r2x);

   r4x=fx(v(i)+h*r3v);    %4ª derivada *1
   r4v=fv(x(i)+h*r3x);

   x(i+1)=x(i)+h*(1/6*r1x+1/3*r2x+1/3*r3x+1/6*r4x);
   v(i+1)=v(i)+h*(1/6*r1v+1/3*r2v+1/3*r3v+1/6*r4v);
end

%% SOLUÇÃO ANALÍTICA

w  = sqrt(k/m);

xt = x(1)*cos(w.*t);
vt = -x(1)*w*sin(w.*t);

%% GRÁFICOS

figure(1)
subplot(2,1,1)
plot(t,v,'b',t,vt,'r')
grid on
title('Comparação Runge-Kutta com solução analitica')
xlabel('tempo (s)')
ylabel('Velocidade (m/s)')
legend('Runge-Kutta','Analitica')

subplot(2,1,2)
plot(t,x,'b',t,xt,'r')
grid on
title('Comparação Runge-Kutta com solução analitica')
xlabel('tempo (s)')
ylabel('Posição (m)')
legend('Runge-Kutta','Analitica')