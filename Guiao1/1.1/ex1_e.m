%%  Ex1.e

clear,clc,close all

%%  CONSTANTES

m=0.150;        %150g => unidades SI: 0.150 kg
g=9.8;      %aceleraçao gravitica
v0=0;       %velocidade inicial
alt=3*2;        %altura de um andar*nºde andares que é largada a pedra
t0=0;       %tempo inicial
tf=2;       %tempo final
h=0.2;      %passo temporal

t=t0:h:tf;  %vetor tempo

%%  MÉTODO DE EULER

v=[];
v(1)=v0;
i=1;

while i<length(t)
    v(i+1)=v(i)+(-g)*h;
    i=i+1;
end

%%  GRÁFICO

figure(1)
plot(t,v)
grid on
title('Estimativa numérica da velocidade')
xlabel('tempo (s)')
ylabel('velocidade (m/s)')