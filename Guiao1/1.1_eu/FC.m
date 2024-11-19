t0 = 0;
tf = 2;
h = 0.25;
t = t0:h:tf;
N = length(t);
v = zeros(N,1);
v(1) = v0

for k=1:N
 v(k+1) = v(k) + 1 ;
end

plot(t,v)