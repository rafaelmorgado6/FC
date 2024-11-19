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
    
    % Runge-Kutta 4ª Ordem; func Anónimas
    f_fPSI = @(VV,PSI,EE) 2*PSI*(VV-EE);
    f_PSI = @(fPSI) fPSI;

    for TOL = [1e-3 1e-6 1e-9]
        
        tic
        for is = 1:Nmax
        
            psi = zeros(1,N); psi(end) = psi0;
            fPsi = zeros(1,N); fPsi(end) = fPsi0;
            V = x;
        
            % Runge-Kutta 4ª Ordem
            for k = N-1:-1:1
                r1fPSI = f_fPSI(V(k+1),psi(k+1),E(is));
                r1PSI = f_PSI(fPsi(k+1));
                
                r2fPSI = f_fPSI(V(k+1)+1/2*h,psi(k+1)+1/2*r1PSI*h,E(is));
                r2PSI = f_PSI(fPsi(k+1)+1/2*r1fPSI*h);
                
                r3fPSI = f_fPSI(V(k+1)+1/2*h,psi(k+1)+1/2*r2PSI*h,E(is));
                r3PSI = f_PSI(fPsi(k+1)+1/2*r2fPSI*h);
                
                r4fPSI = f_fPSI(V(k+1)+h,psi(k+1)+r3PSI*h,E(is));
                r4PSI = f_PSI(fPsi(k+1)+r3fPSI*h);
            
                fPsi(k) = fPsi(k+1)-(1/6*r1fPSI+1/3*r2fPSI+1/3*r3fPSI+1/6*r4fPSI)*h;
                psi(k) = psi(k+1)-(1/6*r1PSI+1/3*r2PSI+1/3*r3PSI+1/6*r4PSI)*h;
            end
            %
            
            PF(is) = psi(1);
            
            
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