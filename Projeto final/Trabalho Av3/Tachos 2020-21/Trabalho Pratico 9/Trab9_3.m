clc 
close all
clear all


% Mudar para encontrar outros valores proprios:
E(1)=3;
E(2)=2;

a=1; % Metade da largura do poço
b=4*a;% Mudar para encontrar outros valores proprios:

% Potencial
V0=10; % Quanto nos estendemos para os lados
h=0.001;
x=-b:h:b;
N=numel(x);

% Shooting and matching
tolera=1e-10;
nmaxE=100; % Numero maximo de iteracoes
x_match=0.5;%desde que seja diferente de zero e convem estar entre -a e a
ind_match=round(1+(x_match+b)/h);
x_left=x(1:ind_match);
x_right=x(ind_match:N);
N_left=numel(x_left);
N_right=numel(x_right);

% Potencial
V_left=zeros(1,N_left);
V_left(x_left<-a)=V0;
V_right=zeros(1,N_right);
V_right(x_right>a)=V0;



% Primeiros dois valores de psi 
psi_extremo=0;
psi_seguinte_left=h/100000; % Tanto faz?
psi_seguinte_right=psi_seguinte_left;
% psi_seguinte_right=-psi_seguinte_left; % Solucoes pares/impares

figure(1)
xlabel('x');ylabel('y')



for iE=1:nmaxE
    
    % Esquerda para a direita
    % Constantes auxiliares para o metodo de Numerov
    g=2*(E(iE)-V_left); 
    aux1=(1+h^2/12*g);
    aux2=2*(1-5*h^2/12*g);
    
    psi_left=zeros(1,N_left);
    psi_left(1)=psi_extremo;
    psi_left(2)=psi_seguinte_left;
          
    for n=2:N_left-1
        psi_left(n+1)=(-aux1(n-1)*psi_left(n-1)+aux2(n)*psi_left(n))/aux1(n+1);
    end
    
    % Direita para a esquerda
    clear g
    % Constantes auxiliares para o metodo de Numerov
    g=2*(E(iE)-V_right); % Neste problema, k ao quadrado  muda
    aux1=(1+h^2/12*g);
    aux2=2*(1-5*h^2/12*g);
    
    psi_right=zeros(1,N_right);
    psi_right(N_right)=psi_extremo;
    psi_right(N_right-1)=psi_seguinte_right;
          
    for n=N_right-1:-1:2
        psi_right(n-1)=(-aux1(n+1)*psi_right(n+1)+aux2(n)*psi_right(n))/aux1(n-1);
    end
    
    plot(x_left,psi_left,x_right,psi_right)
    xlabel('x');ylabel('\psi')
    pause (1.0)
    
    % Acerta os valores
     ratio=psi_left(N_left)/psi_right(1);
     psi_right=psi_right*ratio;
%     psi_left=psi_left/sqrt(abs(ratio));
%     psi_right=psi_right*sqrt(abs(ratio))*sign(ratio);
    
    plot(x_left,psi_left,x_right,psi_right)
    xlabel('x');ylabel('\psi')
    pause(1.0) % Para ir vendo como nos aproximamos da solucao.
      
    
    D_left=(25/12*psi_left(N_left)-4*psi_left(N_left-1)+3*psi_left(N_left-2)...
            -4/3*psi_left(N_left-3)+1/4*psi_left(N_left-4))/h;
    D_right=(-25/12*psi_right(1)+4*psi_right(2)-3*psi_right(3)...
            +4/3*psi_right(4)-1/4*psi_right(5))/h;
   
        
   DLog_left=D_left/psi_left(N_left);
   DLog_right=D_right/psi_right(1);
   result(iE)=(DLog_left-DLog_right)/(DLog_left+DLog_right);
            
           
    if(iE>1)
        m=(result(iE)-result(iE-1))/(E(iE)-E(iE-1));
        if m==0; % Sem isto pode dar erro para tolerancias pequenas
            break
        end    
        E(iE+1)=E(iE)-result(iE)/m;
        
        if (abs(E(iE+1)-E(iE)) < tolera)
            break
        end
        
    end
end

% E1_exato=pi^2/8/a^2
fprintf('Energia = %f Ha\n', E(end))


psi=[psi_left(1:N_left-1), psi_right];

% Normalizacao
c_norm=sqrt(trapz(x,psi.^2));
psi_norm=psi/c_norm;
figure(2);
plot(x,psi_norm)
xlabel('{\it x}');ylabel('\psi_{norm}')
