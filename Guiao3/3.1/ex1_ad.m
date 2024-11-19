%%  Ex1.a

% Equação diferencial principal:
%             Fr=m*a <=>
%             <=> m*a=-k*x <=>
%             <=> m*v'=-k*x <=>
%             <=> m*x'' = -k*x <=>
%             <=> x'' = -k*x/m
% 
%  Dividindo esta EDO em 2 de primeira ordem temos:
%             v' = -k*x/m     e     v=x'


%       0   |
%       1/2 | 1/2
%       ---------------
%           | 0     1

%%  EX1.d
clear,clc,close all

%% CONSTANTES

h=0.01;
tf=10;
t=0:h:tf;

k=16;
m=1;

x=[];
v=[];

x(1)=1;
v(1)=0;

%% MÉTODO DE RUNGE-KUTTA

%Equações das derivadas:  (Só é necessário colocar as equações)
fx = @(t,X,V) V;        %derivada em ordem ao tempo de x
fv = @(t,X,V) -k*X/m;   %derivada em ordem ao tempo de v

for i=1:length(t)-1
                %VER PARTE DE CIMA ( DE CIMA PARA BAIXO)
    r1x= fx(t(i)+0*h,x(i),v(i));    %primeiro valor =0 ou seja,
    r1v= fv(t(i)+0*h,x(i),v(i));    %r1x = fx(t(i) +0*h, x(i) +0*h*r0x, ...)    

                                %segundo valor =1/2 ou seja, h*1/2
    r2x= fx(t(i)+h/2 , x(i)+r1x*h/2 , v(i)+r1v*h/2);
    r2v= fv(t(i)+h/2 , x(i)+r1x*h/2 , v(i)+r1v*h/2);

                    %VER PARTE DE BAIXO DA TABELA(esquerda r1x -> direita
                    %r2x)
    x(i+1)=x(i)+ r2x*h; %x(i) + h* ( 0*r1x + 1*r2x)
    v(i+1)=v(i)+ r2v*h; %v(i) + h* ( 0*r1x + 1*r2x)
end

%% SOLUÇÃO ANALÍTICA

w  = sqrt(k/m);

xt = x(1)*cos(w.*t);
vt = -x(1)*w*sin(w.*t);

%% GRÁFICOS

figure(1)
subplot(2,1,1)
plot(t,v,'b',t,vt,'r')
grid on
title('Comparação Runge-Kutta com solução analitica')
xlabel('tempo (s)')
ylabel('Velocidade (m/s)')
legend('Runge-Kutta','Analitica')

subplot(2,1,2)
plot(t,x,'b',t,xt,'r')
grid on
title('Comparação Runge-Kutta com solução analitica')
xlabel('tempo (s)')
ylabel('Posição (m)')
legend('Runge-Kutta','Analitica')