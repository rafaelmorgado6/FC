clc; clear;  close all,
% Inicialização de variáveis
xmax = 10;
h = 0.001;
x = 0:h:xmax;
N = length(x);

Nmax = 500; TOL = 1e-9;
B = x(1); 
E(1) = 1.7; E(2) = 1.8;

%ciclo iterativo para shooting
jj = 0;
EE = zeros(1,3);

for k = [1.8 3.1 4.0]
    jj = jj + 1;
    E = [k  k+0.1];
    for is = 1:Nmax

        g(2:N) = -2.*(x(2:N) - E(is));  
        C1 = (1+h^2/12.*g);  C2 = (1-5*h^2/12.*g);

        PSI = zeros(1,N);
        PSI(N) = 0; PSI(N-1) = h;

        % Método Numerov regressivo
        for k = N-1:-1:2
            PSI(k-1) = C1(k-1)^(-1)*(2*C2(k)*PSI(k) - C1(k+1)*PSI(k+1));
        end
        %

        PF(is) = PSI(1);

        if is>1
           %shooting 
            declive = (PF(is)-PF(is-1))/(E(is)-E(is-1));

            if declive == 0
                break
            end

            E(is+1) = E(is) + (B-PF(is))/declive;

            if(abs(E(is)-E(is+1)) <= TOL)
                break
            end
            %

        end
    end
 EE(jj) = E(end)
 psiE(jj,:) = PSI;
end


%% Valores próprios c)
figure(2)
h = 0.001;
x = -10:h:10;

ai = airy(x);

plot(x,ai)
hold on
yline(0)
N = length(x);

i = 0;
for k = N:-1:2
    if or(and(ai(k)>0,ai(k-1)<=0),and(ai(k)<=0,ai(k-1)>0))
        i = i+1;
        an(i) = interp1(ai(k-1:k),x(k-1:k),0,'linear');
    end
    if i == 3
        break
    end
end

En = -2^(-1/3)*an

%% funções próprias e)
for n = 1:3
    a= airy((2^(1/3).*(x-En(n))));
    fp(n,:) = a(10001:end);
end

