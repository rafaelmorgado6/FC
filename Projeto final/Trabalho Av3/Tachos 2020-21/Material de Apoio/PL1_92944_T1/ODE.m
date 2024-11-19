
function derivadas = ODE(t,solution,m,K,alpha,mu,w0,F0)
  
    derivadas = zeros(2,1);
    y=solution(1);
    v=solution(2);
    
    derivadas(1)= v;
    derivadas(2)=(-K*y*(1+alpha*y^2) +mu*cos(v)*v + F0*cos(w0*t))/m;
end