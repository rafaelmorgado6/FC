
clear all
close all
clc


%constantes
L=2;
tol = 10^(-13); %tolerância para o shooting
h = 0.1;


E1 = pi^2/(2*L^2); % energia para o nível n=1: objetivo 

x = 0:h:L;
Nx = length(x);
psi = zeros(1,Nx); %função de onda
v=zeros(1,Nx);v(Nx)=1; %dpsi/dx- 1º derivada, arbitrar um valor inicial

% guesses
E(1) = E1 + 1;
E(2) = E1 - 1;

% Condições fronteira
psi(1)=0;
psi(Nx) = 0;
psi(Nx-1) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov

% Goal do shooting
B = 0; % psi_f = 0

k = 1;
comparator = 1;

fx=@(V) V; %primeira derivada
fv=@(E, psi) -2*E*psi; %segunda derivada


while abs(B - comparator) > tol %método shooting
    
    %runge-kutta 4ºordem
    for i=Nx:-1:2  
      r1x=fx(v(i));
      r1v=fv(E(k),psi(i));

      r2x=fx(v(i)+r1v*h/2);
      r2v=fv(E(k)+r1x*h/2,psi(i)+r1x*h/2);

      r3x=fx(v(i)+r2v*h/2);
      r3v=fv(E(k)+r2x*h/2,psi(i)+r2x*h/2);

      r4x=fx(v(i)+r3v*h);
      r4v=fv(E(k)+r3x*h/2,psi(i)+r3x*h/2);

      psi(i-1)=psi(i)+(r1x/6+r2x/3+r3x/3+r4x/6)*h;
      v(i-1)=v(i)+(r1v/6+r2v/3+r3v/3+r4v/6)*h;
    end
    
    comparator = psi(1);
    psi_f(k) = comparator;
    
    if k > 1
        m = (psi_f(k)-psi_f(k-1))/(E(k)-E(k-1));
        E(k+1) = E(k) + (B-psi_f(k))/m;
    end
    k = k + 1;
end


% Normalização
C = trapz(x,abs(psi).^2);
psi_norm = psi/sqrt(C);
plot(x,psi_norm)
xlabel('posição,x')
ylabel('psi')
title('função de onda normalizada')
fprintf('o valor obtido para a energia em n=1, através dos métodos, é: %s .\n',num2str(E(k)));


