
clear all
close all
clc


%constantes
L=2;
tol = 10^(-13); %tolerância para o shooting
h = 0.1;
C=5*pi^2/(2*L^3);
ET=[pi^2/(2*L^2), 2*pi^2/L^2, 9*pi^2/(L^2*2)] % energia para o nível n=1,2,3 teórica 

x = 0:h:L;
Nx = length(x);
psi = zeros(3,Nx); %função de onda %cada linha para cda nivel
g=zeros(3,Nx);


% Condições fronteira
psi(:,1)=0;
psi(:,2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov


for j=1:3 %ciclo para calcular cada valor de energia de cada nivel
    
    E1=ET(j)
    
    % guesses
    E(j,1) = E1 + 2;
    E(j,2) = E1 - 2;
    
    % Goal do shooting
    B = 0; % psi_f = 0
    
    k = 1;
    comparator = 1;

    while abs(B - comparator) > tol %método shooting

        g(j,:)= 2*(E(j,k)-C*x);

        % Método Numerov progressivo
        for i = 2:Nx-1
            psi(j,i+1) = (1+h^2*g(j,i+1)/12)^(-1)*(-(1+h^2*g(j,i-1)/12)*psi(j,i-1) + 2*(1-5*h^2*g(j,i)/12)*psi(j,i));
        end

        comparator = psi(j,end);
        psi_f(j,k) = comparator;

        if k > 1
            m = (psi_f(j,k)-psi_f(j,k-1))/(E(j,k)-E(j,k-1));
            E(j,k+1) = E(j,k) + (B-psi_f(j,k))/m;
        end
        k = k + 1;
        
    end
    
end
    C1 = trapz(x,abs(psi(1,:)).^2);
    C2 = trapz(x,abs(psi(2,:)).^2);
    C3 = trapz(x,abs(psi(3,:)).^2);
    
    psi_norm1 = psi(1,:)/sqrt(C1);
    psi_norm2 = psi(2,:)/sqrt(C2);
    psi_norm3 = psi(3,:)/sqrt(C3);
    
    figure(1)
    plot(x,psi_norm1)
    xlabel('posição,x')
    ylabel('psi')
    title('função de onda normalizada de E1')
    
    figure(2)
    plot(x,psi_norm2)
    xlabel('posição,x')
    ylabel('psi')
    title('função de onda normalizada de E2')
    
   figure(3)
    plot(x,psi_norm3)
    xlabel('posição,x')
    ylabel('psi')
    title('função de onda normalizada de E3')

    

    