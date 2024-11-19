% Fábio Caldas, 80248, P4
% Inês Leite, 98490, P4


clc; clear; close all;

% Inicialização de variáveis
xmax = 10;
h = 0.001;
x = 0:h:xmax;
N = length(x);

Nmax = 500; TOL = 1e-9;
B = x(1); 
E(1) = 1.8; E(2) = 1.7;

%ciclo iterativo para shooting
tic
for is = 1:Nmax
    
    g(2:N) = -2.*(x(2:N) - E(is));  
    C1 = (1+h^2/12.*g);  C2 = (1-5*h^2/12.*g);

    PSI = zeros(1,N);
    PSI(N) = 0; PSI(N-1) = h;
    
    % Método Numerov regressivo
    for k = N-1:-1:2
        PSI(k-1) = C1(k-1)^(-1)*(2*C2(k)*PSI(k) - C1(k+1)*PSI(k+1));
    end
    %
    
    PF(is) = PSI(1);
        
    if is>1
        %shooting 
        declive = (PF(is)-PF(is-1))/(E(is)-E(is-1));
        
        if declive == 0
            break
        end
        
        E(is+1) = E(is) + (B-PF(is))/declive;
        
        if(abs(E(is)-E(is+1)) <= TOL)
            break
        end
        %
        
    end
end
tempo = toc

fprintf('E1 = %f \n', E(end)) %Energia no estado fundamental

% Normalização
Cnorma = sqrt(trapz(x,PSI.^2));
PSI_normal= PSI/Cnorma;

% Gráficos
plot(x,PSI_normal) 
xlabel('x');ylabel('\psi_{norm}')
title("\psi_{norm} em função de x")