%%  Ex2.d

clear,clc,close all

%% CONSTANTES

tf=10;      %TROCAR AQUI O TEMPO FINAL

k=16;
m=1;

x=[];
v=[];

x(1)=1;
v(1)=0;

%% MÉTODO DE RUNGE-KUTTA de 4ª ordem
%PARA SIMPLIFICAR SÓ É NECESSÁRIO COLOCAR AS VARIAVEIS DAS EXPRESSOES, NÃO
%É NECESSÁRIO COLOCAR (t,X,V)

fx = @(V) V;        %derivada em ordem ao tempo de x
fv = @(X) -k*X/m;   %derivada em ordem ao tempo de v

h=[0.01 0.02 0.05 0.2 0.1];

for j=1:length(h)
    t=0:h(j):tf;

    for i=1:length(t)-1 
        r1x=fx(v(i));    %1ª derivada *0
        r1v=fv(x(i));

        r2x=fx(v(i)+h/2*r1v);    %2ª derivada *1/2
        r2v=fv(x(i)+h/2*r1x);
    
        r3x=fx(t(i)+h/2, v(i)+h*0*r1v+h/2*r2v);    %3ª derivada *1/2
        r3v=fv(x(i)+h/2*r2x);

        r4x=fx(v(i)+h*r3v);    %4ª derivada *1
        r4v=fv(x(i)+h*r3x);

        x(i+1)=x(i)+h*(1/6*r1x+1/3*r2x+1/3*r3x+1/6*r4x);
        v(i+1)=v(i)+h*(1/6*r1v+1/3*r2v+1/3*r3v+1/6*r4v);
    end

    erro(j)= abs(x(N) - x_tfin_sa);
end
%% GRAFICO DO ERRO