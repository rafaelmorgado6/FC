clc; close all; clear;

L = 1; %cm
T0 = 0; %�C
TL = 100; %�C
alpha = 0.01; %cm^2/s
u = 0.1; %cm/s

tf = 50; dt = 0.125; t = 0:dt:tf; %s
Nt = length(t);
dx = 0.05; x = 0:dx:L; %cm
Nx = length(x);

T_t0 = 100.*x./L;

T = zeros(Nx,Nt); T(:,1) = T_t0; T(1,:)= T0; T(Nx,:) = TL;

D = alpha*dt/(dx^2); C = u/2*dt/dx;

for n = 1:Nt-1 %ciclo do tempo
    for i = 2: Nx-1 %ciclo do comprimento
        T(i,n+1) = T(i,n) + D*(T(i+1,n)-2*T(i,n)+T(i-1,n))-C*(T(i+1,n)-T(i-1,n));
    end
end

figure(1); plot(x,T(:,1:10:end))
figure(2);mesh(t,x,T)
figure(3)
Tx = 100.*(exp(10.*x./L)-1)./(exp(10)-1);
plot(x,Tx);
