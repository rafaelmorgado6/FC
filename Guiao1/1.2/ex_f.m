%%  Ex2.f

clear,clc,close all

%%  CONSTANTES

Vc0=5;      %Tensão inicial condensador
C=1e-3;     %Condensador
L=0.25;     %Bobina
t0=0;       %tempo em que o S fecha
tf=0.5;     %tempo final
h=0.0001;    %passo temporal

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

%% SOLUÇÃO ANALÍTICA

Q0=Vc0*C;
w=sqrt(1/(L*C));    %frequencia angular das oscilações

Vct = Vc0*cos(w*t);
Qt=Q0*cos(w*t);
It=-Q0*w*sin(w*t); %ou It=-Imax*sin(w*t), em que Imax=-w*Q0

%%  GRÁFICO

figure(1)
subplot(3,1,1)
plot(t,Vc,t,Vct)
grid on
title('Comparação do método de Euler com a solução analitica -> Vc')
xlabel('tempo (s)')
ylabel('Volts (V)')
legend('Método de Euler','Solução Analítica')

subplot(3,1,2)
plot(t(1:5000),I,t(1:5000),It(1:5000))
grid on
title('Comparação do método de Euler com a solução analitica -> I')
xlabel('tempo (s)')
ylabel('Corrente (A)')
legend('Método de Euler','Solução Analítica')

subplot(3,1,3)
plot(t(1:5000),Q,t(1:5000),Qt(1:5000))
grid on
title('Comparação do método de Euler com a solução analitica -> Q')
xlabel('tempo (s)')
ylabel('Carga (C)')
legend('Método de Euler','Solução Analítica')

%Quanto menor o passo temporal mais exato é o metodo de euler
%comparadamente a sua solução analitica do problema em questão. Em termos
%de amplitude esta quanto mais pequeno o passo temporal mais igual a
%amplitude do método de euler será a solução analitica. Quanto ao periodo
%este também será igual com a diminuição do passo temporal