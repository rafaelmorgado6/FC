%%  Ex2.g

clear,clc,close all

%%  GRÁFICO PARA CALCULAR O PERÍODO

Vc0=5;      %Tensão inicial condensador
C=1e-3;     %Condensador
L=0.25;     %Bobina
t0=0;       %tempo em que o S fecha
tf=0.5;     %tempo final
h=0.001;    %passo temporal

a=1/(L*C);  %constante

t=t0:h:tf;  %vetor tempo

DV=[];
Vc=[];

DV(1)=0;
Vc(1)=Vc0;

for i=1:length(t)-1
    DV(i+1)=DV(i)-a*Vc(i)*h;    %METODO DE EULER PARA A DERIVADA DE Vc
    Vc(i+1)=Vc(i)+DV(i)*h;      %METODO DE EULER PARA Vc 
end

figure(1)
plot(t,Vc)
grid on
title('Tensão em funçao do tempo')
xlabel('tempo (s)')
ylabel('Volts (V)')

w=sqrt(1/(L*C));    %frequencia angular das oscilações

%%  CÁLCULO DO PERÍODO PRÁTICO E TEÓRICO

ind=find(islocalmax(Vc)); %islocalmax serve para criar um vetor de 0 e 1 
                          %que apenas existem 1 nos máximos. Depois o find
                          %é para saber o índice em que estes 1 se situam
                      
tt=(t(ind));         %para ver os valores dos tempos correspondentes aos 
                     %indices já achados

for i=2:length(tt)
    Ppra(i-1)=tt(i)-tt(i-1);    %calculo de todos os períodos possiveis
end

Tpratico=mean(Ppra) %PERIODO PRATICO

Tteorico=2*pi/w     %PERIODO TEORICO