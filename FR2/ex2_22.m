%% alinea 2.2
clear all
close all
clc

kmola = 2;
m =1.5;

a = kmola/m;
w =sqrt(a);

T = (2*pi())/w;

alfa_(1) = -0.18;
alfa_(2) = -0.2;
t0 = 0;
tf = 50;
t = t0:0.01:tf;
x0 = 1.9;
v0 = 0;


AN = -1.5;
B = AN;


%tolerancia
reltol = 3e-14;
abstol_1 =1e-13;
abstol_2 =1e-13;

for is=1:2000

    alfa = alfa_(is);
    options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
    [t,solucao] = ode45(@f,[t0 tf],[x0 v0],options,kmola,m,alfa);

    %soluçõesk
    x = solucao(:,1);
    v = solucao(:,2);



imin = 0;

    for k=2:length(t)-1
        if((x(k+1)>x(k)) && (x(k-1) >= x(k)))
            imin = imin+1;
            aux = lagr(t(k-1:k+1),-x(k-1:k+1));
            xmin(imin)=-aux(2);
        end
    end
    
    AN(is)= mean(xmin);
    
    if(is>1)
        declive = (AN(is)-AN(is-1))/(alfa_(is)-alfa_(is-1));
        alfa_(is+1)= alfa_(is) + (B-AN(is))/(declive);
        
        if(abs(B-AN(is))<10^(-4))
            break
        end
    end
end
figure
plot(x)
figure
plot(t,x);