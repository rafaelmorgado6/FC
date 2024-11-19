 %%  Ex1.e

clear,clc,close all

%% EQUAÇÕES DIFERENCIAS PARA APLICAÇÃO DO MÉTODO DE EULER

% Equação diferencial principal:
%             Fr=m*a <=>
%             <=> m*a=-k*x <=>
%             <=> m*v'=-k*x <=>
%             <=> m*x'' = -k*x <=>
%             <=> x'' = -k*x/m
% 
%  Dividindo esta EDO em 2 de primeira ordem temos:
%             v' = -k*x/m     e     v=x'

%% CONSTANTES

x0=1;   
v0=0;
k=1;
m=1;

h=0.01;     %passo temporal (IMPOSTO POR NÓS)
t0=0;

w=sqrt(k/m);    %calculo do tempo final, uma aproximação do tempo final
T=2*pi/w;
tf=T*8;

t=t0:h:tf;      %vetor tempo

v=[];
x=[];
v(1)=v0;
x(1)=x0;

%% MÉTODO DE CRACK-NICOLSON

% AZ=b ;   A=[1 -h/2 ; w^2*h/2 1]  e b=[x(i)+h/2*v(i) ; v(i)-w^2*h/2*x(i)]

A=[1 -h/2; w.^2*(h/2) 1];

for i=1:length(t)-1
    B=[x(i)+(h/2)*v(i); v(i)-w.^2*(h/2)*x(i)];

    Z=linsolve(A,B);
    x(i+1)=Z(1);
    v(i+1)=Z(2);
end

%% SOLUÇÃO ANALÍTICA

%NÃO DÃO A SOLUÇÃO ANALÍTICA COMO É QUE QUEREM QUE COMPARE-MOS

%% ENERGIA MECÂNICA

Em=m*v.^2*0.5+k*x.^2*0.5;   %Em=Ec+Epe= 1/2*m*v^2 + k*x^2/2

%% CALCULO DA AMPLITUDE E DO PERÍODO DE MOVIMENTO

imax=0;

for i=2:length(t)-1
    if ((x(i+1)-x(i)<=0) && (x(i)-x(i-1)>=0)) %condição x(i-1)<=x(i)=>x(i+1)

            imax = imax + 1;
            aux=lagr(t(i-1:i+1),x(i-1:i+1)); %maximos locais
            tmax(imax)=aux(1)
            xmax(imax)=aux(2);
    end
end 

p=polyfit(1:imax,tmax,1);
Tpratico=p(1)              % t = 6.2884
Apratico=mean(xmax)        % amplitude = 1.0126 -> visto que é Euler-Cromer amplitude é aproximadamente 1

Aanalitica=sqrt(x(1).^2+(m*v(1).^2)/k)  % Aanalitica=1
Tanalitica=2*pi/w                       % Tanalitica=6.2832