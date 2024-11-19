clear all
close all
clc


%constantes
L=2;
tol = 10^(-13); %toler�ncia para o shooting
h = 0.05;
b=-0.06152;
dt=5;

ET=pi^2/(2*L^2); % energia para o n�vel n=1: objetivo 

x = 0:h:L;
Nx = length(x);
psi = zeros(1,Nx); %fun��o de onda
g=zeros(1,Nx);

Ne = 100; %numero de itera��es do shooting
E = zeros(1,Ne); %energia 
psi_f = zeros(1,Ne); %fun��o de onda final

C=5*pi^2/(2*L^3);
% Condi��es fronteira
psi(1)=0;
psi(2) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov



    % guesses
    E1=ET
    E(1) = E1 + 2;
    E(2) = E1 - 2;
    % Goal do shooting
    B = 0; % psi_f = 0
    
    k = 1;
    comparator = 1;

    while abs(B - comparator) > tol

        g(1,:)= 2*(E(k)-C*x);


        % M�todo Numerov progressivo
        for i = 2:Nx-1
            psi(i+1) = (1+h^2*g(i+1)/12)^(-1)*(-(1+h^2*g(i-1)/12)*psi(i-1) + 2*(1-5*h^2*g(i)/12)*psi(i));
        end

        comparator = psi(end);
        psi_f(k) = comparator;

        if k > 1
            m = (psi_f(k)-psi_f(k-1))/(E(k)-E(k-1));
            E(k+1) = E(k) + (B-psi_f(k))/m;
        end
        k = k + 1;
        
    end
    

    C1 = trapz(x,abs(psi).^2);
    psi_norm1 = psi/sqrt(C1)
   
   cte=((pi/dt)^(2/3).*(dt.*x/L-E(k)/ET));
   Ai=airy(cte)
   Bi=airy(2,cte)
psi1=Ai+b*Bi;

C = trapz(x,abs(psi1).^2);
psi1_norm=psi1/sqrt(C);


plot(x,psi1_norm,'b',x,psi_norm1,'r')
xlabel('posi��o,x')
ylabel('psi Et')
title('fun��o de onda normalizada')
    













