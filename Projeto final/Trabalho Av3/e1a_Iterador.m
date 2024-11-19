% Fábio Caldas, 80248, P4
% Inês Leite, 98490, P4


clc; clear; close all;

% Inicialização de vetores
it = [];
tempo = [];
dif_h = [];

for h = [0.01 0.001 0.00001]

    % Inicialização de variáveis
    xmax = 10;
    x = 0:h:xmax;
    N = length(x);

    psi0 = 0;
    fPsi0 = -h;
    
    Nmax = 500;
    B = x(1);
    E(1) = 1.8; E(2) = 1.7;

    for TOL = [1e-3 1e-6 1e-9]
        
        tic
        for is = 1:Nmax
        
            g(2:N) = -2.*(x(2:N) - E(is));  
            C1 = (1+h^2/12.*g);  C2 = (1-5*h^2/12.*g);
        
            PSI = zeros(1,N);
            PSI(N) = 0; PSI(N-1) = -h;
            
            % Método Numerov regressivo
            for k = N-1:-1:2
                PSI(k-1) = C1(k-1)^(-1)*(2*C2(k)*PSI(k) - C1(k+1)*PSI(k+1));
            end
            %
            
            PF(is) = PSI(1);
            
            
            if is>1
               % shooting 
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
        % Armazenamento dos tempos e números de iterações
        it = [it is];
        tempo = [tempo toc];
    end
    % Armazenamento dos vetores ocupados anteriormente e posterior desocupação
    dif_h = [dif_h; it; tempo];
    it = [];
    tempo = [];
end