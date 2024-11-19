%% Autores:
% André da Silva Santos, 92933, PL9
% Fábio Pinto Caldas, 80248, PL9

%% e)

clc; clear all; close all;

xmax = 10;
h = 0.001;
x = 0:h:xmax; % tentativa-erro
N = numel(x);
nmaxE = 100;

    % Shooting
tol = 1e-12;

% Energia aproximada 
nivel = input(['Escolha o valor de n\n']);
fprintf('\n')
E_aprox = (((3*pi)/sqrt(2))*(nivel-(1/4)))^(2/3);

E(1) = 0.99*E_aprox;     % guess(0)
E(2) = 1.01*E_aprox;     % guess(1)


    % Numerov
psi = zeros(1,N);
psi(N) = 0;
psi(N-1) = h/1000; 


for iE=1:nmaxE
    
    g(2:N) = 2*E(iE)-4*x(2:N);
    
    % Constantes auxiliares para realização do metodo de Numerov
    % Regressivo:
    aux1 = (1+((h^2)/12)*g);
    aux2 = 2*(1-(5*(h^2)/12)*g);
    
    % Método de Numerov, aplicado regressivamente
    for k = (N-1):-1:3
        psi(k-1) = (-aux1(k+1)*psi(k+1) + aux2(k)*psi(k))/aux1(k-1);
    end
    psi(1) = interp1(x(2:5),psi(2:5),0,'spline');  
    
    plot(x,psi)
    xlabel('x'); ylabel('psi');
    
    
    % Shooting
    psi_f(iE) = psi(1);    
    
    if(iE>1)
        m = (psi_f(iE) - psi_f(iE-1))/(E(iE) - E(iE-1));
        if m==0; 
            break
        end    
        E(iE+1) = E(iE) - psi_f(iE)/m;
       
        if abs(psi(1)-0) < tol
            break
        end
    end
end

fprintf('E(end) = %f Ha\n', E(end))
fprintf('E_aprox = %f Ha\n', E_aprox)
fprintf('E_aprox/E(end) = %f \n', E_aprox/E(end))

% Normalização
c_norm = trapz(x,psi.^2);
psi_norm = psi./sqrt(c_norm);
figure;
plot(x,psi_norm,'y*')
xlabel('{\it x}');ylabel('{\it X}_{norm}')

% Probabilidade
Energia = E(iE);
xo = Energia/2;
prob = trapz(x(round(xo/h):end),psi_norm(round(xo/h):end).^2);
fprintf('Probabilidade de encontrar a particula numa região classicamente proibida = %f \n', prob)