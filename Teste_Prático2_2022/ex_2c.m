%%  Ex2c

clear,clc,close all

%%  CONSTANTES

g=9.8;
m=0.5;
DD=0:0.5/9:0.5;  %10 valores

h=0.1;
t=0:h:20;
N=length(t);

Vx=length(N);
Vy=length(N);
x=length(N);
y=length(N);

y(1)=20;
x(1)=0;

B=50;   %objetivo do problema
tol=10e-4;
v0(1)=40;
v0(2)=41;

%%  METODO DE SHOOTING
j=1;
for D=0: 0.5/9 : 0.5

for i=1:100

    %velocidade inicial muda com o tempo logo tem de estar dentro do ciclo
    Vx(1)=v0(i)*cosd(10);
    Vy(1)=v0(i)*sind(10);

    for k=1:N-1     %APLICAÇÃO DO EULER
        Vx(k+1)=Vx(k)-D/m*Vx(k)*h;
        Vy(k+1)=Vy(k)+(-g-D/m*Vy(k))*h;
        
        x(k+1)=x(k) + Vx(k)*h;
        y(k+1)=y(k) + Vy(k)*h;

        if(y(k+1)<0)
            y=y(1:k+1);     %CORTAR VETORES PARA INTERP
            x =x(1:k+1);
            Vy=Vy(1:k+1);
            Vx=Vx(1:k+1);
           break;
        end
    end

    xsolo = interp1(y(end-1:end),x(end-1:end),0,'linear');
    xend(i)=xsolo;

    if(i>1)     
        %METODO DA SECANTE
        m=(xend(i)-xend(i-1))/(v0(i)-v0(i-1));    %formula do declive normal
                                                %result=x    e    guess=v0

        v0(i+1)=v0(i)+(B-xend(i))/m;      %guess(3) = guess(2) + (B-result(2))/m                                
    
        if(abs(xend(i)-B)<tol)   %condição de saida
          break;
        end
    end
end

    yD(j,:)=y;
    xD(j,:)=x;
    v0D(j)=v0(end);
    j=j+1;
end
%% GRÁFICOS

figure(1)
for l=1:10
    plot(xD(l,:),yD(l,:))
    hold on
end
title('y em função de x')
xlabel('x')
ylabel('y')
grid on

figure(2)
plot(DD,v0D)
xlabel('D')
ylabel('v0')
grid on