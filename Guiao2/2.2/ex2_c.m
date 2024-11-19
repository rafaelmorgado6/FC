%%  Ex2.c

clear,clc,close all

%% EQUAÇÕES DIFERENCIAS PARA APLICAÇÃO DO MÉTODO DE EULER

% Equação diferencial principal:
%             Fr=m*a <=>
%             <=> m*ar=-G*ms*mm*rr/r^3 <=>
%             <=> vr'=-G*rr/r^3 <=>
%             <=> rr'' = -G*rr/r^3 <=>
%             <=> (xî+yî)'' = -G*(xî+yî)/r^3 <=>
%             <=> xî'' = -G/r^3*xî   e     yî'' = -G/r^3*yî
% 
%  Dividindo esta EDO em 2 de primeira ordem temos:
%             vr' = -G*rr/r^3     e     vr=rr'

%% CONSTANTES

G=4*pi^2;       %Produto Gms(constante gravitacional)
h=0.0001;       %Passo temporal (ano)
t=0:h:1;        %vetor tempo

%USAR PLANO XY PARA REPRESENTAR A FORÇA, logo teremos posição em x e em y e
%velocidades de x e velocidades de y !!SOL É A ORIGEM DO REFERENCIAL!!

x=[];
y=[];
vx=[];
vy=[];
r=[];       %raio da trajetoria em uma posição exata -> para calcular a força

x(1)=0.47; 
y(1)=0;
vx(1)=0;
vy(1)=8.2;
r(1)=sqrt(x(1)^2+y(1)^2);   %raio no momento 

%% MÉTODO DE EULER-CROMER

for i=1:length(t)-1
    vx(i+1)=vx(i)-G/r(i)^3*x(i)*h;   %velocidade em x
    x(i+1)=x(i)+vx(i+1)*h;      %posição x

    vy(i+1)=vy(i)-G/r(i)^3*y(i)*h;   %velocidade em x
    y(i+1)=y(i)+vy(i+1)*h;      %posição x

    r(i+1)=sqrt(x(i+1)^2+y(i+1)^2); %raio da posição seguinte
    
end

%% CALCULO DA AREA VARRIDA

theta=[];
Area=[];
Area(1)=0;
theta(1)=0;
for i=1:length(t)-1
    theta(i+1)=mod(atan2(y(i+1),x(i+1)),2*pi); %obter coordenada polar theta
    Area(i+1)=r(i+1)^2*(theta(i+1)-theta(i))/2;

    if ((y(i+1)<y(1)) && (y(i)>y(1)))   %condiçao para meia orbita
        break
    end
end

%% GRÁFICOS
                %%  TA MAL CARALHOOOOOOO
figure(1)
plot(t(2:length(Area)),Area(2:length(Area)))
grid on
title('Variação da Area varrida em funçao do tempo')
xlabel('Tempo (anos)')
ylabel('Variação Area (AU^2)')