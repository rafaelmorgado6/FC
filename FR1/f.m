function derivadas = f(t,solucao,m,K)
% A função f retorna as derivadas dos valores do tempo e da solução da
% função f
derivadas=zeros(2,1); % alocar as derivadas
x=solucao(1);
v=solucao(2);
% O vetor soulção tem os valores de x e v para o tempo t em que a função é
% chamada pela rotina de ode45
% Condição inicial é dada pelo programa principal
derivadas(1)=v;
derivadas(2)=-K/m*x;
end