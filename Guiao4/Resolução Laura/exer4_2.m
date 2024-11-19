close all;
clc;
clear all;

u=10^-3;
L=1;
T=10^3;

x0=0;
n=5;
xf=L;
h=0.01;

x=x0:h:xf;
N=length(x);
n=N-2;

A=eye(n,n);         %Matriz identidade
A=-2*A;
%METODO DAS DIFERENÇAS FINITAS
for k=1:n
    A(k,k+1) = 1;       %vai fazer o triangulo sup
    A(k+1,k) = 1;       %vai fazer o triangulo inf
    %toda a matriz vai ficar a '1' menos a diagonal que ficaria -2
end

[vec,val] = eigs(A,3,'sm')      %vai dar 2 matrizes: vec e val
                                %onde val sao os valores da frequencia
                                %fundamental
                                %vec são os ys - variação da corda

vals = diag(val);               %vai me dar os valores proprios para 
                                % calcular as frequencias

w=sqrt(-vals*T/u)/h             %calculo das frequencias

m1=vec(:,1)';                 
m2=vec(:,2)';
m3=vec(:,3)';

subplot(3,1,1)
plot(x(1:end-1),m1,"r")
title('1º Modo de Vibração')

subplot(3,1,2)
plot(x(1:end-1),m2,"r")
title('2º Modo de Vibração')

subplot(3,1,3)
plot(x(1:end-1),m3,"r")
title('3º Modo de Vibração')
