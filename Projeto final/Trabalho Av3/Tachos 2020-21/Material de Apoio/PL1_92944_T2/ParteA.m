%Ema Fadiga PL1-92944

clc
clear all
close all
format long

%Alínea b

%Constantes
K=1;
M=1;
b=0.2;
B=1.3;% É o valor de result que se pretende

%Condições iniciais
y0=0;
v0=-1.5;

%Vetor t
tf=60;
t=0:0.05:tf;
Nt=numel(t);

%Valores para guess
guess=[1.4,1.6];

tol=1e-10;
it_max=50;%número de iterações máximo para que o ciclo não se repita
options = odeset('RelTol',3E-14,'AbsTol',[1E-13 1E-13]);

for i=1:it_max
    mu=guess(i);
    %aplicação da ode45
    [t,solution] = ode45(@funode,t,[y0 v0],options,K,M,b,mu);
    y=solution(:,1);
    Ny=numel(y);
        
    %Cáluclo dos máximos (Amplitudes positivas)
    imax=0;
    for k=500:Ny-1
        if (y(k)>y(k-1) && y(k)>y(k+1))
            imax=imax+1;
            aux=lagr(t(k-1:k+1),y(k-1:k+1));
            ymax(imax)=aux(2);
            tmax(imax)=aux(1);
        end
    end
    %Amplitude
    result(i)=mean(ymax);
   
    %Método da secante
    if(i>1)
        m=(result(i)-result(i-1))/(guess(i)-guess(i-1));
        guess(i+1)=guess(i)+(B-result(i))/m;
              
        if (abs(result(i)-B)<tol)
            mu_obt=guess(i);
            break
        end
    end    
end
Result_obt=result(end)%Verifica-se que é igual a B    
u=mu_obt%Valor obtido para mu

%Alinea c)
[t,solution] = ode45(@funode,t,[y0 v0],options,K,M,b,mu);
y=solution(:,1);
v=solution(:,2);

imax=0;
for j=500:Nt-1
    if (y(j)>y(j-1) && y(j)>y(j+1))
            imax=imax+1;
            aux=lagr(t(j-1:j+1),y(j-1:j+1));
            ymax(imax)=aux(2);
            tmax(imax)=aux(1);
            
            if imax>1
                T(imax-1)=tmax(imax)-tmax(imax-1);
            end       
    end
end
Periodo=mean(T)

N1=round(Nt/3);

figure (1)
plot(t,y)
title('y(m)')
xlabel('t(s)')
ylabel('y(m)')

figure (2)
plot(y(N1:2*N1),v(N1:2*N1))
title('Espaço de fases')
xlabel('y(m)')
ylabel('v(m/s)')
