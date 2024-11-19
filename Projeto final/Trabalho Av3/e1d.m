% Fábio Caldas, 80248, P4
% Inês Leite, 98490, P4


clc; clear;  close all;

xmax = 10; h = 0.001;
x = 0:h:xmax;
N = length(x);

Nmax = 500; TOL = 1e-9;
B = 0; 
E(1) = 1.7; E(2) = 1.8;

count = 0; EE = zeros(1,3);

%ciclo iterativo para shooting
for j = [1.8 3.1 4.0]
    count = count + 1;
    E = [j  j+0.1];
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
    EE(count) = E(end);
    Cn = sqrt(trapz(x,PSI.^2));
    PSIn_pratico(count,:) = PSI./Cn;
end

% Gráficos
plot(x,PSIn_pratico)
xlabel('x');ylabel('\psi_{norm}')
title("\psi_{norm} em função de x")
grid on; legend('E_1','E_2','E_3');