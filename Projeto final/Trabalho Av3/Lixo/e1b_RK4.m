clc; clear all; close all;

% Constantes
E = 1.8;

x0 = 0; xf = 10;
h = 0.001;

psi0 = 0;
fpsi0 = -h;

% Vetor das Posições
x = x0:h:xf;
N = length(x);

% Inicialização de vetores
psi = zeros(1,N); psi(end) = psi0;
fPsi = zeros(1,N); fPsi(end) = fPsi0;
V = x;

% Runge-Kutta 4ª Ordem
f_fPSI = @(VV,PSI) 2*PSI*(VV-E);
f_PSI = @(fPSI) fPSI;

for k = N:-1:2
    r1fPSI = f_fPSI(V(k),psi(k));
    r1PSI = f_PSI(fPsi(k));
    
    r2fPSI = f_fPSI(V(k)+1/2*h,psi(k)+1/2*r1PSI);
    r2PSI = f_PSI(fPsi(k)+1/2*r1fPSI);
    
    r3fPSI = f_fPSI(V(k)+1/2*h,psi(k)+1/2*r2PSI);
    r3PSI = f_PSI(fPsi(k)+1/2*r2fPSI);
    
    r4fPSI = f_fPSI(V(k)+h,psi(k)+r3PSI);
    r4PSI = f_PSI(fPsi(k)+r3fPSI);

    fPsi(k-1) = fPsi(k)+(1/6*r1fPSI+1/3*r2fPSI+1/3*r3fPSI+1/6*r4fPSI)*h;
    psi(k-1) = psi(k)+(1/6*r1PSI+1/3*r2PSI+1/3*r3PSI+1/6*r4PSI)*h;
end

plot(x,psi);