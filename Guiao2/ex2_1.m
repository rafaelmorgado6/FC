% a)
clc;
clear all;
close all;


m = 1;
kmola = 1;
h = 0.1;
t0 = 0;
tf = 50;
t = t0:h:tf;
n = length(t);
x = zeros(1,n);
vx = zeros(1,n);
Dvx = zeros(1,n);
vx(1) = 0;
x(1) = 1;

for i=1:(n-1)
  vx(i+1) = vx(i) + (-kmola/m)*x(i)*h;
  x(i+1) = x(i) + vx(i)*h;
end

for i=1:n
    Et(i) = (1/2)*m*vx(i).^2+1/2*kmola*x(i).^2;
end


plot(t,Et)
title('Método de Euler')

% b) 
clc;
clear all;
close all;


m = 1;
kmola = 1;
h = 0.1;
t0 = 0;
tf = 50;
t = t0:h:tf;
n = length(t);
x = zeros(1,n);
vx = zeros(1,n);
Dvx = zeros(1,n);
vx(1) = 0;
x(1) = 1;

for i=1:(n-1)
  vx(i+1) = vx(i) + (-kmola/m)*x(i)*h;
  x(i+1) = x(i) + vx(i+1)*h;
end

for i=1:n
    Et(i) = (1/2)*m*vx(i).^2+1/2*kmola*x(i).^2;
end


plot(t,Et)
title('Método de Euler-Cromer')

% c)

clc;
clear all;
close all;


m = 1;
kmola = 1;
h = 0.1;
t0 = 0;
tf = 50;
t = t0:h:tf;
n = length(t);
x = zeros(1,n);
vx = zeros(1,n);
Dvx = zeros(1,n);
vx(1) = 0;
x(1) = 1;

for i=1:(n-1)
  vx(i+1) = (vx(i)-x(i)*h)/(1+h^2);
  x(i+1) = x(i) + vx(i+1)*h;
end

for i=1:n
    Et(i) = (1/2)*m*vx(i).^2+1/2*kmola*x(i).^2;
end


plot(t,Et)
title('Método de Euler Implícito')

% d)
   
clc;
clear all;
close all;


m = 1;
kmola = 1;
h = 0.01;
t0 = 0;
tf = 50;
t = t0:h:tf;
n = length(t);
x = zeros(1,n);
vx = zeros(1,n);
Dvx = zeros(1,n);
vx(1) = 0;
x(1) = 1;
A = [1 h/2; -h/2 1];

for i = 1:n-1
    b = [vx(i) - h/2*x(i); x(i)+vx(i)*h/2];
    sol = linsolve(A,b);
    vx(i+1) = sol(1);
    x(i+1) = sol(2);
end

for i=1:n
    Et(i) = (1/2)*m*vx(i).^2+1/2*kmola*x(i).^2;
end


plot(t,Et)
title('Método de Crank–Nicolson')


