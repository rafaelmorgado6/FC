%%  Ex3.a

clear,clc,close all

%% CONSTANTES

m=1;
k=1;
a=-0.1;

h=0.01;    %h random nossa decisão 
w=sqrt(k/m);    %calculo do tempo final, uma aproximação do tempo final
T=2*pi/w;
tf=T*8;

t=0:h:tf;      %vetor tempo

x=[];
v=[];

x(1)=1;
v(1)=1;

%% METODO DE EULER-CROMER

for i=1:length(t)-1
    v(i+1)=v(i)-k/m*(x(i)+2*a*x(i)^3)*h;
    x(i+1)=x(i)+v(i+1)*h;
end

%% ENERGIA MECANICA

Em=m*v.^2*0.5+k*x.^2*0.5;   %Em=Ec+Epe= 1/2*m*v^2 + k*x^2/2

%% CALCULO DA AMPLITUDE E DO PERIODO DE OSCILAÇÃO

counter=1;
for i=1:length(t)-1
    if (i > 1) & ((x(i) > x(i+1)) & (x(i)>x(i-1)))
        x_max(counter)=x(i);
        t_max(counter)=t(i);
        ind(counter) = i;
        counter = counter +1;
    end
end

for i = 1:length(ind)
    aux = lagr(t(ind(i)-1:ind(i)+1),x(ind(i)-1:ind(i)+1));
    t_lagr(i) = aux(1);
    x_lagr(i) = aux(2);
end

A = mean(x_lagr);
T = t_max(2)-t_max(1);

%% GRÁFICOS

figure(1)
subplot(3,1,1)
plot(t,v)
grid on
title('Método de Euler-Cromer velocidade em funçao do tempo')
xlabel('tempo (s)')
ylabel('Velocidade (m/s)')

subplot(3,1,2)
plot(t,x)
grid on
title('Método de Euler-Cromer posição em funçao do tempo')
xlabel('tempo (s)')
ylabel('Posição (m)')

subplot(3,1,3)
plot(x,v)
grid on
title('Método de Euler-Cromer velocidade em funçao da posição')
xlabel('Posição (m)')
ylabel('Velocidade (m/s)')

figure(2)
plot(t,Em)
grid on
title('Euler-Cromer => Energia mecânica em função do tempo')
xlabel('tempo (s)')
ylabel('Energia mecânica (J)')