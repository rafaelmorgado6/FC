
%% i
clear all
close all
clc

U=10^-3; % Kg/m
L=1; % m
T=10^3; % N
x0=0; 
n=1; 
xf=L;
dy0=0; fdy0=0.01;   % não sabemos, assumir que é sempre assim
w=(n*pi/L)*sqrt(T/U);

h=0.0001; 
x=x0:h:xf;
N=length(x);
dy=zeros(1,N);
fdy=zeros(1,N);
dy(1)=dy0;
fdy(1)=fdy0;

C=(w.^2)*U/T;   %dy'
%Metodo euler
for kk=1:N-1
    fdy(kk+1)=fdy(kk)-C*dy(kk)*h;
    dy(kk+1)=dy(kk)+fdy(kk+1)*h;
end
   
figure(1)
plot(x,dy,"r")

%Trata-se sim de um modo fundamental pois começa em '0' e acaba em '0'

%% ii)
clear all
close all
clc

U=10^-3; % Kg/m
L=1; % m
T=10^3; % N
x0=0;
n=1; 
xf=L;
dy0=0; fdy0=0.01;
%w=(n*pi/L)*sqrt(T/U);
w=2000;

h=0.0001; 
x=x0:h:xf;
N=length(x);
dy=zeros(1,N);
fdy=zeros(1,N);
dy(1)=dy0;
fdy(1)=fdy0;

C=(w.^2)*U/T;   %dy'    
%Metodo euler
for kk=1:N-1
    fdy(kk+1)=fdy(kk)-C*dy(kk)*h;
    dy(kk+1)=dy(kk)+fdy(kk+1)*h;
end
   
figure(1)
plot(x,dy,"r")

%neste caso nao se considera uma frequencia fundamental

%% METODO DE SHOOTING
clear,clc,close

U=10^-3; % Kg/m
L=1; % m
T=10^3; % N
x0=0; 
n=5; 
xf=L;
dy0=0; fdy0=0.01;
w=(n*pi/L)*sqrt(T/U);

h=0.0001; 
x=x0:h:xf;
N=length(x);
dy=zeros(1,N);
fdy=zeros(1,N);
dy(1)=dy0;
fdy(1)=fdy0;
w(1)=2000;
w(2)=3500;

tol=1e-10;
B=0; % valor que queremos como valor final
for is=1:100                %100 valores que se vai calcular
    C=(w(is).^2)*U/T;       %vai ser sempre diferente 

    for kk=1:N-1            % Euler
        fdy(kk+1)=fdy(kk)-C*dy(kk)*h;
        dy(kk+1)=dy(kk)+fdy(kk+1)*h;
    end

    dyf(is)=dy(end);  %vai sempre buscar o ultimo valor obtido de y
                      %para depois fazer melhores aproximaçoes

    if (is>1)     %so vai ser verdade na segunda iteração para 
                  %obter sempre novos valores

        if abs(w(is))<tol   %requisito
            break
        end

        %METODO DAS SECANTES
        m=(dyf(is)-dyf(is-1))/(w(is)-w(is-1));
        %b=dyf(is)-m*w(is);
        w(is+1)=w(is)+(B-dyf(is)/m);
    end
    if (abs(w(is+1)-w(is)))<tol     %igual à controlo: se a diferença
                                    %entre as duas amostras for mais
                                    %pequena que a tolerancia, sai
        break
    end
end

w(is)
figure(1)
plot(x,dy,"r")