%%  Ex1c

clear,clc,close all

%%  CONSTANTES

a=0.01;
u=0.1;
L=1;

tf=50;
dt=0.125;
dx1=0.1;
dx2=0.01;

%vetores
x1=0:dx1:L;
x2=0:dx2:L;
t=0:dt:tf;

Nt=length(t);
Nx1=length(x1);
T1=zeros(Nx1,Nt);
Nx2=length(x2);
T2=zeros(Nx2,Nt);

T1(:,1)=100*x1/L; %valores iniciais
T1(1,:)=0;       %valores das faces
T1(Nx1,:)=100;

T2(:,1)=100*x2/L; %valores iniciais
T2(1,:)=0;       %valores das faces
T2(Nx2,:)=100;

D1=a*dt/dx1^2;
C1=u*dt/(2*dx1);

D2=a*dt/dx2^2;
C2=u*dt/(2*dx2);

%%  METODO DE EULER

for n=1:Nt-1
    for i=2:Nx1-1
    
        T1(i,n+1)=T1(i,n) +D1*(T1(i-1,n)-2*T1(i,n)+T1(i+1,n)) - C1*(T1(i+1,n)-T1(i-1,n));

    end
end

for n=1:Nt-1
    for i=2:Nx2-1
    
        T2(i,n+1)=T2(i,n) +D2*(T2(i-1,n)-2*T2(i,n)+T2(i+1,n)) - C2*(T2(i+1,n)-T2(i-1,n));

    end
end

%%  CRITERIO DE ESTABILIDADE

nx1= 2*a/u;
nt1= dx1^2/(2*a);
disp([num2str(dx1), ' <= ', num2str(nx1),' e ' ,num2str(dt), ' <= ',num2str(nt1)])
disp('Logo dx=0.1 cm é estável')

nt2= dx2^2/(2*a);
disp([num2str(dx1), ' <= ', num2str(nx1),' e ' ,num2str(dt), ' <= ',num2str(nt2)])
disp('Logo dx=0.01 cm não é estável')


