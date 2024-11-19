%%  Ex1.d

clear,clc,close all

%%  CONSTANTES

m=0.150;    %150g => unidades SI: 0.150 kg
g=9.8;      %aceleraçao gravitica
v0=0;       %velocidade inicial
alt=3*2;    %altura de um andar*nºde andares que é largada a pedra
t0=0;       %tempo inicial
tf=2;       %tempo final

%%  MÉTODO DE EULER

v=[];       %vetor velocidade
v(1)=v0;    %velocidade inicial do vetor = 0

h=[10e-3, 10e-4, 10e-5, 10e-6, 10e-7];

for j=1:length(h)

    t=t0:h(j):tf;  %vetor tempo

    for i=1:length(t)-1
        v(i+1)=v(i)+(-g)*h(j);     %metodo de euler    
    end

    vz=v0-g*t;
    
    vmin=min(v);      %ver o valor maximo de Vc
    vzmin=vz(end);
    ERROGLOBAL(j)=abs(vzmin-vmin);  %Calculo do erro global para cada h  
end
%%  GRÁFICO

figure(1)
plot(log(h),log(ERROGLOBAL),'*')
grid on
title('Erro Global em função do passo temporal')
xlabel('log(h)')
ylabel('log(ERROGLOBAL)')