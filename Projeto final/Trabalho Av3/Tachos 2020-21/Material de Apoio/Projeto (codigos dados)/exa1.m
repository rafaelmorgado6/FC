
clear all
close all
clc


%constantes
L=2;
tol = 10^(-13); %toler�ncia para o shooting
h = 0.1;

E1 = pi^2/(2*L^2) % energia para o n�vel n=1: objetivo 

x = 0:h:L;
Nx = length(x);
psi = zeros(1,Nx); %fun��o de onda


% guesses
E(1) = E1 + 1; %primeiras guesses tendo em conta o valor te�rico
E(2) = E1 - 1;

% Condi��es fronteira para usar o m�todo de Numerov Regressivo
psi(1)=0;
psi(Nx) = 0;
psi(Nx-1) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov

% Goal do shooting
B = 0; % psi_f = 0

k = 1;
comparator = 1; 

while abs(B - comparator) > tol %m�todo de shooting
    
    g= 2*E(k);
    const=(h^2/12)*g;
    
    % M�todo Numerov regressivo
    for i = Nx-1:-1:2
        psi(i-1) = (1+const)^-1 *(2*psi(i)*(1-5*const) - (1+const)*psi(i+1));
    end
    
    comparator = psi(1);
    psi_f(k) = comparator; %valor de psi usado no m�todo da secante
    
    if k > 1
        m = (psi_f(k)-psi_f(k-1))/(E(k)-E(k-1));
        E(k+1) = E(k) + (B-psi_f(k))/m;
    end
    k = k + 1;
end

% Normaliza��o
C = trapz(x,abs(psi).^2);
psi_norm = psi/sqrt(C);

plot(x,psi_norm)

xlabel('posi��o,x')
ylabel('psi')
title('fun��o de onda normalizada')
fprintf('o valor obtido para a energia em n=1, atrav�s dos m�todos, �: %s .\n',num2str(E(k)));



