clc;
clear all;
close all; 

t0 = 0;
tf = 1.6;
h = 0.2;
t = t0:h:tf;
N = length(t);
v0 = 0;
m = 0.150;
g = 9.8;
v = zeros(N,1);
z = zeros(N,1);
v(1) = 0;
z(1) = 6;

% analitico
for k=2:N-1
    za = z(1) - 0.5*g*t.^2;
    if za(k+1)<0
        break
    end
end

% numérico
 for k=1:N-1
    v(k+1) = v(k) + (-g)*h ;
 end

for k=1:N-1
    z(k+1) = z(k) +v(k)*h ;
    if z(k+1)<0
        break
    end
end

plot(t(1:k+1),z(1:k+1),t(1:k+1),za(1:k+1))


tsolo = interp1(z(k:k+1),t(k:k+1),0,'linear')

vsolo = v(1) - g*tsolo
