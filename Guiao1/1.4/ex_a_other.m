%%  MÉTODO DE EULER -> CALCULO DO ERRO GLOBAL

DV=[];
Vc=[];
Q=[];
I=[];
ERROGLOBAL=[];

DV(1)=0;
Vc(1)=Vc0;

h=[10e-1, 10e-2, 10e-3, 10e-4];

for j=1:length(h)

    t=t0:h(j):tf;  %vetor tempo

    for i=1:length(t)-1
        DV(i+1)=DV(i)-a*Vc(i)*h(j);    %METODO DE EULER PARA A DERIVADA DE Vc
        Vc(i+1)=Vc(i)+DV(i)*h(j);      %METODO DE EULER PARA Vc 

        I(i)=DV(i)*C;   %CALCULO DA CORRENTE 
        Q(i)=Vc(i)*C;    %CALCULO DA CARGA
    end
    
    Vcmax= max(Vc);      %ver o valor maximo de Vc
    ERROGLOBAL(j)=abs(Vc0-Vcmax);  %Calculo do erro global para cada h  
end

figure(1)
plot(log(h),log(ERROGLOBAL),'b')
grid on
title('Erro Global em função do passo temporal')
xlabel('log(h)')
ylabel('log(ERROGLOBAL)')
