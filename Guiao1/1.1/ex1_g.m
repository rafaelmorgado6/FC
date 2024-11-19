%%  Ex1.g

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

%%  GRÁFICO

figure(1)
plot(t(1:length(z)-1),z(1:length(z)-1)) %para não mostrar valores negativos
grid on
title('Posição em Função do tempo')
xlabel('tempo (s)')
ylabel('Posição (m)')

%%  CALCULO DO INSTANTE EM QUE A PEDRA CAI NO CHÃO E SUA VELOCIDADE

tsolo=interp1(z(i:i+1),t(i:i+1),0,'linear') %estimativa numerica do instante 
                                            %em que a pedra cai no chão
                                            %PARAMETROS INTERP1:(queremos o
                                            %tempo)
                                            %1-> ALTURA
                                            %2-> TEMPO
                                            %3-> POSIÇÃO DO ACONTECIMENTO=0
                                            %4-> LINEAR                              
 
vsolo=v0-g*tsolo    %velocidade que atinge o solo