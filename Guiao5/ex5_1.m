%%  Ex5.1

close,clc,clear all

%% Constantes 

L=50;   %cm

k=0.93;
c=0.094;
p=8.9;

tf=500;
dx=0.5;         %variar para ver estabilidade
dt=0.1;         %variar para ver estabilidade

%vetores
t=0:dt:tf;
x=0:dx:L;

Nt=length(t);
Nx=length(x);
T=zeros(Nx,Nt);

T(:,1)=100; %tempo = inicial (1) -> temperatura = 100
T(1,:)=0;  %colocar o valor inicial da extremidade em 0
T(Nx,:)=0;  %colocar a outra extremidade = 0

D=k/(c*p);  %Constante criada so para simplicar a escrita das equações

%%  Método de Euler

for n=1:Nt-1
    for i=2:Nx-1
    
        T(i,n+1)=T(i,n) +D *1/dx^2 *(T(i-1,n)-2*T(i,n)+T(i+1,n)) * dt; %dt é o passo temporal

    end
end

%%  Gráficos

figure(1)
contourf(x,t,T')
xlabel('x')
ylabel('t')
zlabel('ºC')

figure(2)
mesh(t,x,T)
xlabel('t')
ylabel('x')
zlabel('ºC')

%%  Critério de estabilidade

%D=k/(c*p);
n= D * dt/dx^2;  %se for menor ou igual que 1/2 é estável
disp(['É estável pois: ', num2str(n) , ' <= 1/2']);

%%  Ver tempo que demoram os pontos a L/4 a diminuir até 50ºC

L4 = L/4;   %ver quais os pontos que temos de ver

index_L4 = ((Nx-1)/4 + 1);  %descobrir qual o indice de x a que está L/4    
T_xL4(:,1)= T(index_L4, :); %T(index_L4,:) é o vetor da temperatura em x=L/4
                            %T_xL4(:,1) = colocar os valores da temperatura
                            %num novo vetor so para se fazer a interpolação
                            %mais facilmente

figure(3)
plot(t, T(index_L4,:))
xlabel('t(s)')
ylabel('T(ºC)')
grid on
title('Temperatura a x=L/4')

t_50=interp1(T_xL4(T_xL4<90), t(T_xL4<90),50);  %achar o tempo (é preciso as condiçoes 
                                                %para nao haver repetiçoes
                                                %de valores e não causar
                                                %erros)
disp(['tempo para T=50ºC a x=L/4: ', num2str(t_50)])    %so para aparecer no terminal o valor