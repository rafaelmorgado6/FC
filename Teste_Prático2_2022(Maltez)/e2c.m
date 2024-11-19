clc; clear; close all;


ang = 10; %º
ym = 20;  %m
m = 0.5; %Km
g= 9.8; %m/s^-2
DD = 0:0.01:0.5; %N/ms^(-1)

h = 0.1;t =0:h:20; %s
N = length(t);


B = 50;
TOL = 1e-4;
vi = [41 40];
j = 1;
v0 = zeros(1,length(0:0.05:0.5));

for D = 0:0.01:0.5
    vi = [41 40];
for i = 1:100
    y = zeros(1,N); y(1) = ym;
    vy = zeros(1,N); vy(1) = sind(ang)*vi(i);
    x = zeros(1,N); x(1) = 0;
    vx = zeros(1,N); vx(1) = cosd(ang)*vi(i);
    B = 50;
    
    for k = 1 : N-1
        vy(k+1) = vy(k) - ((D/m)*vy(k) + g )*h;
        y(k+1) = y(k) + vy(k)*h;
        
        vx(k+1) = vx(k) - (D/m)*vx(k)*h;
        x(k+1) = x(k) + vx(k)*h;
        
        if y(k+1) < 0
            y=y(1:k+1);x = x(1:k+1);
            vy=vy(1:k+1);vx=vx(1:k+1);
            break
        end
    end
    
    x_solo = interp1(y(end-1:end),x(end-1:end),0,'linear');
    
    xF(i) = x_solo;
    
    if abs(x_solo - B)<TOL
        break
    end
    if i > 1
        declive = (xF(i) - xF(i-1))/(vi(i)-vi(i-1));
        vi(i+1) = vi(i) + (B - xF(i))/declive;
    end
    
end
    v0(j) = vi(end);
    j = j+1;
end

plot(DD,v0)
title("V0(D)")
xlabel("D / N/ms^(-1)")
ylabel("V0 / m/s")