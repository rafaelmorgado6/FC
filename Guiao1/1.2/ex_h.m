%%  Ex2.h

clear,clc,close all

%%  GRÁFICO PARA CALCULAR A ENERGIA 

Vc0=5;      %Tensão inicial condensador
C=1e-3;     %Condensador
L=0.25;     %Bobina
t0=0;       %tempo em que o S fecha
tf=0.5;     %tempo final
h=0.0001;    %passo temporal

a=1/(L*C);  %constante

t=t0:h:tf;  %vetor tempo

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

figure(1)
plot(t,Vc)
grid on
title('Tensão em funçao do tempo')
xlabel('tempo (s)')
ylabel('Volts (V)')

Q0=Vc0*C;           %Carga inicial

%%  CALCULO DA ENERGIA

Uc=1/2*Q.^2/C;      %Energia Eletrica armazenada no condensador
Um=1/2*L*I.^2;      %Energia Magnetica armazenada na bobina

Utpratico=mean(Uc+Um)

Utteorico=1/2*Q0^2/C

%Como o Utpratico=Utteorico a energia total não se conserva com o método de
%Euler