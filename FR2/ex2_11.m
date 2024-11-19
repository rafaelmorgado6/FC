clc;clear all;close all;
beta=18;
h=0.01;
x=[0:h:7];
Nshooting = 2000;
result=[];

v1=1.9;
v2=2.0;
tolera=1e-4;

guess = [v1 v2];
B=0;

for i = 1:Nshooting
    V = guess(i);
    T=zeros(1,length(x));
    T(1)=1e-4;
    Tlinha=zeros(1,length(x));
    Tlinha(1)=1e-4*V;
for n = 1:length(x)-1
        Tlinha(n+1) = Tlinha(n)+(V*Tlinha(n)-beta*exp(-1/T(n))*(1+Tlinha(n)-V*T(n)))*h;
        T(n+1) = T(n)+Tlinha(n+1)*h;
end
plot(x,T)
pause(0.5)

result = [result Tlinha(length(x))];

if i>1
        mSec = (result(i) - result(i-1))/ (guess(i)-guess(i-1));
        b = guess(i) + (B - result(i))/mSec;
        guess = [guess b];  
    end
    
    if abs(guess(i+1)-guess(i)) < tolera
        V=guess(i);
        break
    end
end
figure();   
plot(x,T);
xlabel('x')
ylabel('T')