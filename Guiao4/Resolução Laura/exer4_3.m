clear,clc,close

%a) Expressão geral??

%b) a temp à final da resistencia=20 T'(0)=0


T_amb=20;
lambda=0.1;
Q=2.1e6;
r0=0;
R=0.001;
h=1e-6;

r=r0:h:R;
N=length(r);
T=zeros(N,1);

A=-2*eye(N);
A(N,N)=1;

for i=2:N-1
    A(i,i+1)=1+h/(2*r(i));
    A(i,i-1)=1-h/(2*r(i));
end
A(1,2)=2;
b=-h^2*Q/lambda*ones(N,1);
b(N)=T_amb; % correct last element

T=linsolve(A,b);
plot(r,T)