clear all, close all, clc

% Constantes físicas
L = 50; %cm
K = 0.93; %cal/(s cm ºC)
C = 0.094; %cal/(g ºC)
RHO = 8.9; %g/cm^3
ALPHA = K/(RHO*C);

% Pré-alocação de vetores
t0 = 0; %s
tf = 500; %s
x0 = 0; %s
xf = L;

dt = 0.1; %s
dx = 0.5; %cm
t = t0:dt:tf;
x = x0:dx:xf; 

Nt = length(t);
Nx = length(x);

T = zeros(Nx, Nt);

% Condições iniciais
T(2:Nx-1,1) = 100;

% Condições fronteira
T(1,:) = 0;
T(Nx,:) = 0;

BETA = K*dt/(C*RHO*dx^2);

A = (2/BETA + 2)*eye(Nx-2);
b = zeros(Nx-2,1);

% Método de Crank-Nicolson
for  i= 1:Nx-2
    if (i > 1)
        A(i,i-1) = -1;
    end
    if (i + 1 < Nx-2)
        A(i,i+1) = -1;
    end
end

DVB = 2/BETA-2;

[L, U, P] = lu(A);
y = zeros(Nx-2,1);
for j = 1: Nt -1
    b = T(1:Nx-2,j)+DVB*T(2:Nx-1,j)+T(3:Nx,j);
    b(1) = b(1) + T(1,j+1);
    b(Nx-2) = b(Nx-2) + T(Nx,j+1);
    y = L\b;
    T(2:Nx-1,j+1) = U\y;
end

figure(1)
mesh(t,x,T)
figure(2)
contourf(x,t,T')
colorbar