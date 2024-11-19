%clc 
close all
clear all

% Para n = 3, l = 1:
rmax=50; n_esperado=3; l=1; E(1)=-0.06; E(2)=-0.07; 

h=0.001;
r=0:h:rmax;
%r(1)=1e-100;
N=numel(r);

% Shooting
tolera=1e-7;
nmaxE=100; % Numero maximo de iteracoes

% Primeiros dois valores de u 
uN=0;
uNmenos1=h/1000; % Tanto faz?
% uN=2*r(N)*exp(-r(N));
% uNmenos1=2*r(N-1)*exp(-r(N-1));

figure(1)
xlabel('r');ylabel('y')

for iE=1:nmaxE
    
    g(2:N)=2*(E(iE)+1./r(2:N))-l*(l+1)./r(2:N).^2; % Neste problema, g  muda
    % Constantes auxiliares para o metodo de Numerov
    aux1=(1+h^2/12.*g);
    aux2=2*(1-5*h^2/12.*g);

    
    u=zeros(1,N);
    u(N)=uN;
    u(N-1)=uNmenos1;
          
    for k=N-1:-1:3
        u(k-1)=(-aux1(k+1)*u(k+1)+aux2(k)*u(k))/aux1(k-1);
    end
    u(1)=interp1(r(2:5),u(2:5),0,'spline');
    
    % u(1:4)
    
    plot(r,u)
    xlabel('r');ylabel('u')
    pause(1.0) % Para ir vendo como nos aproximamos da solucao.
    
    u_f(iE)=u(1);
    
    
    if(iE>1)
        m=(u_f(iE)-u_f(iE-1))/(E(iE)-E(iE-1));
        if m==0; % Sem isto pode dar erro para tolerancias pequenas
            break
        end    
        E(iE+1)=E(iE)-u_f(iE)/m;
        
        %if (abs(E(iE+1)-E(iE)) < tolera)
        if abs(u(1)-0) < tolera
            break
        end
        
    end
end

E1_exato=-0.5/n_esperado^2;
fprintf('Energia = %f Ha\n', E(end))
fprintf('E1_exato/Energia = %f \n', E1_exato/E(end))

% Calcula funcao de onda exceto em r=0;
R(2:N)=u(2:N)./r(2:N);
% Interpola para obter R(1)
R(1)=interp1(r(2:5),R(2:5),0,'spline');

% Normalizacao 
c_norm=sqrt(trapz(r,r.^2.*R.^2));
R_norm=R/c_norm;
figure(2);
plot(r,R_norm,'y*')
xlabel('{\it r}');ylabel('{\it R}_{norm}')
axis([0 rmax -Inf Inf])
hold on

% Plot para orbital 3p
R_exato=-4*sqrt(2)/(27*sqrt(3))*(1-r/6).*r.*exp(-r/3);
plot(r,R_exato,'k','LineWidth',2)

