function derivadas = f(t,solucao,k,m,alfa)
derivadas = zeros(2,1);

x = solucao(1);
v = solucao(2);

derivadas(1) =v;
derivadas(2) =-k/m*x - (3/2)*(k/m)*alfa*x^2;

end