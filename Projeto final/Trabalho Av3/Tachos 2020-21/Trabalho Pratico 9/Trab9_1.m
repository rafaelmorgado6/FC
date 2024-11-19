clc 
close all
clear all

a=1.; % Metade da largura do poço
h=0.001;
x=-a:h:a;
N=numel(x);

% Shooting
tolera=1e-12;
nmaxE=100; % Numero maximo de iteracoes
E(1)=1.5;       % Mudar os valores seguintes
E(2)=1.4;       % para obter outra solucaoo

% Primeiros dois valores de psi 
psi1=0;
psi2=h; % Tanto faz

figure(1)
xlabel('x');ylabel('\psi')

for iE=1:nmaxE
    
    g=2*E(iE); % Neste problema, g nao varia com x
    % Constantes auxiliares para o metodo de Numerov
    aux1=(1+((h^2)/12)*g);
    aux2=2*(1-(5*(h^2)/12)*g);
    
    psi=zeros(1,N);
    psi(1)=psi1;
    psi(2)=psi2;
          
    for n=2:N-1 %porque ja sabemos que psi(1)=0 e psi(N)=0: condicoes fronteira
        psi(n+1)=(-aux1*psi(n-1)+aux2*psi(n))/aux1;
    end         
    
    plot(x,psi)
    xlabel('x');ylabel('\psi')
%     pause(1.0) % Para ir vendo como nos aproximamos da solucao.

   
    % Shooting
    psi_f(iE)=psi(end);
    
    if(iE>1)
        m=(psi_f(iE)-psi_f(iE-1))/(E(iE)-E(iE-1));
        if m==0; % Sem isto pode dar erro para tolerancias pequenas
            break
        end    
        E(iE+1)=E(iE)-psi_f(iE)/m;
        
        if (abs(E(iE+1)-E(iE)) < tolera)
            break
        end
        
    end
end

E1_exato=pi^2/8/a^2;
fprintf('Energia = %f Ha\n', E(end))
fprintf('Energia/E1_exato = %f \n', E(end)/E1_exato)

% Normalizacao
c_norm=trapz(x,psi.^2);
psi_norm=psi/sqrt(c_norm);
figure(2);
plot(x,psi_norm)
xlabel('x');ylabel('\psi_{norm}')
