%%  EX1.d
clear,clc,close all

%% CONSTANTES

h=0.01;
tf=10;
t=0:h:tf;

k=16;
m=1;

x=[];
v=[];

x(1)=1;
v(1)=0;

%% MÉTODO DE RUNGE-KUTTA

%SIMPLIFICAÇÃO -> REMOVER TODAS AS DERIVADAS NÃO NECESSÁRIAS
%Equações das derivadas:  (Só é necessário colocar as equações)
fx = @(V) V;        %derivada em ordem ao tempo de x
fv = @(X) -k*X/m;   %derivada em ordem ao tempo de v

for i=1:length(t)-1
    r1x= fx(v(i));
    r1v= fv(x(i));

    r2x= fx( v(i)+r1v*h/2 );
    r2v= fv( x(i)+r1x*h/2 );

    x(i+1)=x(i)+ r2x*h;
    v(i+1)=v(i)+ r2v*h;
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