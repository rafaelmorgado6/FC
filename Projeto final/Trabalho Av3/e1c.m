% Fábio Caldas, 80248, P4
% Inês Leite, 98490, P4


clc; clear; close all;

h = 0.001; x = -10:h:10;
N = length(x);

ai = airy(x);

counter = 0;
for k = N:-1:2
    if ai(k)*ai(k-1)<=0
        counter = counter + 1;
        % Interpolação
        an(counter) = interp1(ai(k-1:k),x(k-1:k),0,'linear');
        xk(counter) = x(k);         % Índice de x(k)
        xk_1(counter) =x(k-1);      % Índice de x(k-1)
    end
    if counter == 3
        break
    end
end
% Energias
En = -2^(-1/3)*an;

% Gráficos
plot(x,ai);hold on
yline(0);text(an,[0 0 0],'\bullet')
title("Ai(x)")
xlabel("x"); ylabel("Ai")
grid on;