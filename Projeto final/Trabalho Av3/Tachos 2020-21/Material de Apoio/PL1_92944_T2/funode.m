%Ema Fadiga PL1-92944
function derivative = funode(t,solution,K,M,b,mu)

y=solution(1);
v=solution(2);

derivative = zeros(2,1);
derivative(1)=v;
derivative(2)=(mu*(1-v^2)*v -K*(y+b*y^3))/M;