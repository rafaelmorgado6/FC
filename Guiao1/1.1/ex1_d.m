%%  Ex1.d

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

v=[];       %vetor velocidade
v(1)=v0;    %velocidade inicial do vetor = 0

for i=1:length(t)-1
    v(i+1)=v(i)+(-g)*h;     %metodo de euler    
end

%%  GRÁFICO

figure(1)
plot(t,v)
grid on
title('Estimativa numérica da velocidade')
xlabel('tempo (s)')
ylabel('velocidade (m/s)')