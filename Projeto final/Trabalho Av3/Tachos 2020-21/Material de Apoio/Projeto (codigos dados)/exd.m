
clear all
close all
clc


%constantes
L=2;
tol = 10^(-13); %tolerância para o shooting
h = 0.1;


ET = pi^2/(2*L^2) % energia para o nível n=1: objetivo 
C=5*pi^2/(2*L^3)


x = 0:h:L;
Nx = length(x);
psi = zeros(2,Nx); %função de onda
g=zeros(1,Nx);

Ne = 100; %numero de iterações do shooting
E = zeros(2,Ne); %energia 
psi_f = zeros(2,Ne); %função de onda final

% guesses
E(:,1) = ET +1;
E(:,2) = ET - 1;

% Condições fronteira
psi(:,1)=0;

psi(:,2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov

% Goal do shooting
B = 0; % psi_f = 0

k = 1;
comparator = 1;

while abs(B - comparator) > tol
    
    g1= 2*E(1,k);
    const=(h^2/12)*g1;
    
    % Método Numerov regressivo
    for i = 2:Nx-1
            psi(1,i+1) = (1+const)^(-1)*(-(1+const)*psi(1,i-1) + 2*(1-5*const)*psi(1,i));
        end

        comparator = psi(1,end);
        psi_f(1,k) = comparator;
    
    if k > 1
        m = (psi_f(1,k)-psi_f(1,k-1))/(E(1,k)-E(1,k-1));
        E(1,k+1) = E(1,k) + (B-psi_f(1,k))/m;
    end
    k = k + 1;
end
E(1,k)
% Normalização
CN = trapz(x,abs(psi(1,:)).^2);
psi_norm = psi(1,:)/sqrt(CN);










    % guesses

E(2,1) = ET + 2;
E(2,2) = ET- 2;
    % Goal do shooting
B = 0; % psi_f = 0
    
k1 = 1;
comparator1 = 1;

while abs(B - comparator1) > tol

    g(1,:)= 2*(E(2,k1)-C*x);


        % Método Numerov progressivo
    for i = 2:Nx-1
        psi(2,i+1) = (1+h^2*g(i+1)/12)^(-1)*(-(1+h^2*g(i-1)/12)*psi(2,i-1) + 2*(1-5*h^2*g(i)/12)*psi(2,i));
    end

    comparator1 = psi(2,end);
    psi_f(2,k1) = comparator1;

    if k1 > 1
       m = (psi_f(2,k1)-psi_f(2,k1-1))/(E(2,k1)-E(2,k1-1));
       E(2,k1+1) = E(2,k1) + (B-psi_f(2,k1))/m;
   end
   k1 = k1 + 1;
        
end
E(2,k)

C1 = trapz(x,abs(psi(2,:)).^2);
psi_norm1 = psi(2,:)/sqrt(C1);

plot(x,psi_norm,'b',x,psi_norm1,'r')