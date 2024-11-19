%%  Ex1b

clear,clc,close all

%%  CONSTANTES

a=0.01;
u=0.1;
L=1;

tf=50;
dt=0.125;
dx=0.05;

%vetores
x=0:dx:L;
t=0:dt:tf;

Nt=length(t);
Nx=length(x);
T=zeros(Nx,Nt);

T(:,1)=100*x/L; %valores iniciais
T(1,:)=0;       %valores das faces
T(Nx,:)=100;

D=a*dt/dx^2;
C=u*dt/(2*dx);

%%  METODO DE EULER

for n=1:Nt-1
    for i=2:Nx-1
    
        T(i,n+1)=T(i,n) +D*(T(i-1,n)-2*T(i,n)+T(i+1,n)) - C*(T(i+1,n)-T(i-1,n));

    end
end

%%  Gráficos Usar o mesh ou contourf

figure(1)
mesh(t(1:10:end),x,T(:,1:10:end))
xlabel('x')
ylabel('t')
zlabel('T')


%perfil final de T em função de x   TA ERRADO ESTA PARTE 
figure(2)
mesh(t,x,T)
xlabel('x')
ylabel('t')
zlabel('T')


%%  ESTADO ESTACIONÁRIO

Test=100*(exp(10*x/L)-1)/(exp(10)-1);

figure(3)
plot(x,Test)
xlabel('x')
ylabel('T estacionário')