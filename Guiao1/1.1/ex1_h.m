%%  Ex1.h

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
z=[];
v(1)=v0;
z(1)=alt;

for i=1:length(t)-1
    v(i+1)=v(i)+(-g)*h;     %metodo de euler velocidade
    z(i+1)=z(i)+v(i)*h;

    if z(i+1)<0         %posição não pode ir para valores negativos
        break;
    end
end

%%  SOLUÇÃO ANALÍTICA

z0=alt;
zt=z0-1/2*g*t.^2;

%%  GRÁFICO

figure(1)
plot(t(1:length(z)-1),z(1:length(z)-1),'b',t(1:length(z)-1),zt(1:length(z)-1),'r') 
grid on
title('Posição em Função do tempo')
xlabel('tempo (s)')
ylabel('Posição (m)')
legend('Método de Euler','Solução Analítica')

%%  CONCLUSOES

%O comportamento não é identico ao calculado em e) uma vez que neste caso o
%método de Euler não é tão certo pois este usa valores anteriores para
%calcular o valor esperado seguinte, deste modo a solução analítica é mais
%real do que a solução pelo método de Euler