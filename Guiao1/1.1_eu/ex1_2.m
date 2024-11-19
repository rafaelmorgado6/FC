clc;
clear all;
close all; 

h = 0.001;
tf = 0.5;
t = 0:h:tf;
N = length(t);
w = 1/sqrt(0.25*10^-3);
C = 1*10^-3;
L = 0.25;
V = zeros(1,N);
DV = zeros(1,N);
Q = zeros(1,N);
DV(1) = 0;
V(1) = 5;
Q(1) = 0;

% numérico
for k=1:N-1
    DV(k+1) = DV(k) + ((-w^2)*V(k))*h;
    V(k+1) = V(k) + DV(k)*h;
end


Qn = C * V;
In = C * DV;
Vcn = Qn * C;

for k=1:N-1
    Q(k) = V(k) * C;
end

% analitico
for k=1:N
    Vc(k) = V(1)*cos(w*t(k));
    Q(k) = Q(1)*cos(w*t(k));
    I(k) = -Q(1)*w*sin(w*t(k));


end

figure(1)
subplot(2,2,1)
plot(t,Qn,t,Q)
title('Carga')
subplot(2,2,2)
plot(t,In,t,I)
title('Corrent')
subplot(2,2,3)
plot(t,Vcn,t,Vc)
title('Tensão')

% g)
p = find(islocalmax(V)>0);
TT = t(p);
nI = length(TT);

for i=2nI
    j = i-1;
    Texp(j) = TT(i)-TT(i-1);
end

Periodo = mean(Texp)
PeriodoTeorico = 2*pi/w

% h)
Uc = 0.5.*(Q.^2)./C;
Um = 0.5.*L.*I.^2;
Ut = Uc + Um;

figure(2)
plot(t,Ut)
