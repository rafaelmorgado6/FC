%% Oscilador nÃ£o harmÃ³nico
clear; clear all; close all;
massa=1.5; K=2; 
x0=1.9; v0=0; 
t0=0; tf=50; t=0;

  
guess=zeros(1,2);
alpha1=-0.1; guess(1)=alpha1; %valor "sorteado"
alpha2=-0.2; guess(2)=alpha2; %valor "sorteado"
tolera = 1e-4; %tolerÃ¢ncia
Nshooting = 1e3; %NÂº tentativas
B = -1.5; %Amplitude = -1.5

for n=1:Nshooting
    alpha = guess(n);
    
    reltol = 3e-14; abstol_1 = 1e-13; abstol_2 = 1e-13;
    options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
    [t,sol] = ode45(@f,[t0 tf],[x0 v0],options, K, massa, alpha);

    pks=0; locs=0; ind=0; aux=0; tmin=0; xmin=0;
    [pks,locs]=f(-sol(:,1));


for i=1:length(locs)
        ind=locs(i);
        aux=lagr(t(ind-1:ind+1),-sol(ind-1:ind+1,1));
        tmin(i)=aux(1);
        xmin(i)=aux(2);
    end 
    Result(n) = -mean(xmin) ; %'-' porque findpeaks sÃ³ usa positivos
    
    plot(t,sol(:,1));
    pause(0.5);
    
    if(n>1) %SÃ³ nÃ£o corre na 1Âª iteraÃ§Ã£o
        %MÃ©todo da secante
        m = (Result(n) - Result(n-1))/(guess(n) - guess(n-1)); %declive
        guess(n+1) = guess(n) + (B - Result(n))/m;
        if (abs(guess(n+1)-guess(n)) < tolera) %condiÃ§Ã£o de paragem
           break;
        end
    end
    
end

figure('Name','Posição do oscilador')
plot(t,sol(:,1)), title('Oscilador'), xlabel('Tempo /s'), ylabel('Posição x/m')
grid on;