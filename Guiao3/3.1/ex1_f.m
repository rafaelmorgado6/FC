%%  EX1.f

clear,clc,close all

%% CONSTANTES

h=0.01;
tf=10;      %TROCAR AQUI O TEMPO FINAL
t=0:h:tf;

k=16;
m=1;

x=[];
v=[];

x(1)=1;
v(1)=0;

%% MÉTODO DE RUNGE-KUTTA

%Equações das derivadas:  (Só é necessário colocar as equações)
fx = @(t,X,V) V;        %derivada em ordem ao tempo de x
fv = @(t,X,V) -k*X/m;   %derivada em ordem ao tempo de v

for i=1:length(t)-1
    r1x= fx(t(i),x(i),v(i));
    r1v= fv(t(i),x(i),v(i));

    r2x= fx(t(i)+h/2 , x(i)+r1x*h/2 , v(i)+r1v*h/2);
    r2v= fv(t(i)+h/2 , x(i)+r1x*h/2 , v(i)+r1v*h/2);

    x(i+1)=x(i)+ r2x*h;
    v(i+1)=v(i)+ r2v*h;
end

%% METODO DE EULER

vt=[];
xt=[];
vt(1)=0;
xt(1)=1;

for i=1:length(t)-1
    vt(i+1)=vt(i)-k*xt(i)/m*h;
    xt(i+1)=xt(i)+v(i)*h;
end

%% CALCULO DA ENERGIA MECANICA

EmR= m*v.^2*1/2+k*x.^2*1/2;   %Runge-kutta
EmE= m*vt.^2*1/2+k*xt.^2*1/2;   %Euler

%% GRÁFICOS

figure(1)
subplot(2,1,1)
plot(t,v,'b',t,vt,'r')
grid on
title('Comparação Runge-Kutta com método de euler')
xlabel('tempo (s)')
ylabel('Velocidade (m/s)')
legend('Runge-Kutta','Euler')

subplot(2,1,2)
plot(t,x,'b',t,xt,'r')
grid on
title('Comparação Runge-Kutta com método de euler')
xlabel('tempo (s)')
ylabel('Posição (m)')
legend('Runge-Kutta','Euler')

figure(2)
plot(t,EmR,'b',t,EmE,'r')
title('Comparação Energia mecanica')
xlabel('tempo (s)')
ylabel('Energia mecanica (J)')
legend('Runge-Kutta','Euler')