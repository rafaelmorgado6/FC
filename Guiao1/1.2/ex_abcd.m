%%  Ex2.a

% Pela aplicação das leis de kirchhoff podemos analisar que:
%             Vc+Vl=0 <=>
%             <=> Q/C + L*I' = 0 <=>
%             <=> I' = -Q/(C*L) <=>
%             <=> Q'' = -Q/(C*L) <=>
%             <=> C*Vc'' = -C*Vc/(C*L) <=>
%             <=> Vc'' = -Vc/(C*L)   => Equação diferencial

%%  Ex2.b
% Sendo a=1/(L*C) então:
%             Vc'' = -Vc/(C*L)   => Equação diferencial obtida em a)
%               <=> Vc'' = -a*Vc

%%  Ex2.c

% Aplicando o algoritmo de Euler a equação Vc'' = -a*Vc temos que:
%            DV=Vc'=Q'/C=I/C
% 
%                 DV'= -a*Vc    e   DV = Vc'

%%  Ex2.d

clear,clc,close all

%%  CONSTANTES

Vc0=5;      %Tensão inicial condensador
C=1e-3;     %Condensador
L=0.25;     %Bobina
t0=0;       %tempo em que o S fecha
tf=0.5;     %tempo final
h=0.001;    %passo temporal

a=1/(L*C);  %constante

t=t0:h:tf;  %vetor tempo

%%  MÉTODO DE EULER

DV=[];
Vc=[];
Q=[];
I=[];

DV(1)=0;
Vc(1)=Vc0;

for i=1:length(t)-1
    DV(i+1)=DV(i)-a*Vc(i)*h;    %METODO DE EULER PARA A DERIVADA DE Vc
    Vc(i+1)=Vc(i)+DV(i)*h;      %METODO DE EULER PARA Vc 

    I(i)=DV(i)*C;   %CALCULO DA CORRENTE 
    Q(i)=Vc(i)*C;    %CALCULO DA CARGA
end

%%  GRÁFICO

figure(1)
subplot(3,1,1)
plot(t,Vc)
grid on
title('Tensão em funçao do tempo')
xlabel('tempo (s)')
ylabel('Volts (V)')

subplot(3,1,2)
plot(t(1:500),I)
grid on
title('Corrente em funçao do tempo')
xlabel('tempo (s)')
ylabel('Corrente (A)')

subplot(3,1,3)
plot(t(1:500),Q)
grid on
title('Carga em funçao do tempo')
xlabel('tempo (s)')
ylabel('Carga (C)')

%Os gráficos da Tensão e da Corrente estão em fase 180º e o grafico da
%Tensão com a Carga está em fase 0º

