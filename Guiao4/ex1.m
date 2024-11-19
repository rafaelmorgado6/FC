%%                                      Trabalho Prático 4
close all;clear all;clc;
%% Exercício 4_1
% definir constantes
u=10^-3;
L=1;
T=10^3;

%solução analítica
n=1;
wk=((n*pi)/L)*sqrt(T/u)

%variáveis
w(1)=2000;
w(2)=3500;
tol=10^(-12);
dy(1)=2*10^-2;
y(1)=0;
h=0.001;

%vetor x
x= 0:10^-3:L;

for i=1:100
    
    %definir funções
    fy=@(dy) dy;
    fdy=@(y) -((w(i).^2*u)/T)*y;
  
    % método RK de 4ª ordem
    for k=1:numel(x)-1
 r1v = fdy(y(k));
        r1y = fy(dy(k));
        r2v = fdy(y(k) + r1y*(h/2));
        r2y = fy(dy(k) + r1v*(h/2));
        r3v = fdy(y(k) + r2y*(h/2));
        r3y = fy(dy(k) + r2v*(h/2));
        r4v = fdy(y(k) + r3y*h);
        r4y = fy(dy(k) + r3v*h);
    
        dy(k+1) = dy(k) + (1/6)*(r1v+2*r2v+2*r3v+r4v)*h;
        y(k+1) = y(k) + (1/6)*(r1y+2*r2y+2*r3y+r4y)*h;
   
    end
    
    yf(i)=y(end);
    
    % método da secante
    if(i>1) 
        m = (yf(i) - yf(i-1))/(w(i) - w(i-1));
        w(i+1) = w(i)-yf(i)/m;
    end
    
   if(abs(w(i+1)-w(i))<=tol)
       wf = w(i+1)
       break;       
   end
   
end

figure(1)
plot(x,y)