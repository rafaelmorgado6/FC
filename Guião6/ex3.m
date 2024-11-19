%%                                  Trabalho Prático 6
clear all;close all;clc;
% Exercício 3
% vetores x e z

dx=0.05;
dz=0.05;
x=-50:dx:50;
z=0:dx:4;
Nx=length(x);
Nz=length(z);

% matriz q
q=zeros(Nz,Nx);
q(1,:)=sech(x);

% vetor k -> Fourier (frequencia)
dk=2*pi/ (Nx*dx);
k=-Nx/2*dk: dk: (Nx/2-1)* dk;

% preparar ode45
abstol=ones(1,Nx);
abstol=1e-9.*abstol;
options= odeset('RelTol',1e-9,'AbsTol',abstol);

% calcular fft e multiplicar por exp -> cond inicial para ode45
qtexp0= fftshift( fft( q(1,:))).* exp( (i.*(k.^2).*z(1)) / 2);

% chamar ode45
[z, qtexp]=ode45(@nonlinear,z,qtexp0,options,Nx,k); % g

% fourier inversa
for i= 1:Nz
    q(i,:)=ifft( ifftshift( qtexp(i,:).* exp( -(1i.*(k.^2).*z(i)) / 2)));
end

% fazer gráfico
figure(1)
mesh(x,z,abs(q).^2)
