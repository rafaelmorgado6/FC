%Ema Fadiga PL1-92944
clc
clear all
close all

%Constantes
M=1;
K=1;
b=0.2;
mu=1.879830725302758;

%Condições iniciais
y0=0;
v0=-1.5;

%Vetor tempo
h=0.1;
Nt=2^10;
t=0:h:(Nt-1)*h;

%Aplicação da ode45
options = odeset('RelTol',3E-14,'AbsTol',[1E-13 1E-13]);
[t,solution] = ode45(@funode,t,[y0 v0],options,K,M,b,mu);

y=solution(:,1);
v=solution(:,2);

%Alinea e)

%Frequência angular de e)
dw=2*pi/(Nt*h);
w=-Nt/2*dw:dw:(Nt/2-1)*dw;

%Transformada de Fourier da posição para e)
Y=fft(y);
Y=fftshift(Y);
DS=h*abs(Y).^2;

figure(1)
plot(w,DS,'b-')
title('Densidade espetral')
xlabel('w')
ylabel('DE')

%Alinea f)

%Frequência angular de f)
dw2=dw/4;
w2=-Nt/2*dw2:dw2:(Nt/2-1)*dw2;
h2=2*pi/(Nt*dw2);
t2=0:h2:(Nt-1)*h2;

[t2,solution2] = ode45(@funode,t2,[y0 v0],options,K,M,b,mu);
y2=solution2(:,1);

%Transformada de Fourier da posição para f

Y2=fft(y2);
Y2=fftshift(Y2);
DS2=h2*abs(Y2).^2;

figure(2)
plot(w2,DS2,'r-')
title('Densidade espetral')
xlabel('w2')
ylabel('DE2')

% % PARA COMPARAÇÃO- Gráficos da Densidade espetral
% figure ()
% subplot(2,1,1);plot(w,DE,'r');ylabel('DE Alínea e)')
% title('Densidade espetral')
% subplot(2,1,2);plot(w2,DE2,'b');ylabel('DE2 Alínea f)');xlabel('w')
    
%DE para valores inferiores a 10
DEinf10=DS(DS<10);
winf10=w(DS<10);

figure(3)
plot(winf10,DEinf10,'b-')
title('Densidade espetral < 10')
xlabel('w')
ylabel('DE')

% DE2inf10=DS(DS2<10);
% w2inf10=w(DS2<10);
% 
% % PARA COMPARAÇÃO
% figure(5)
% plot(w2inf10,DE2inf10,'r-')
% title('Densidade espetral < 10')
% xlabel('w2')
% ylabel('DE2')

% PARA COMPARAÇÃO- Gráficos da Densidade espetral para valores <10
% figure ()
% subplot(2,1,1);plot(winf10,DEinf10,'r');ylabel('DE')
% title('Densidade espetral < 10')
% subplot(2,1,2);plot(w2inf10,DE2inf10,'b');ylabel('DE2');xlabel('w')


%Alinea h

%Aceleração dada pela Transformada de Fourier (é igual para os dois casos)
A=Y.*(i.*transpose(w)).^2;
A=ifftshift(A);
atf=ifft(A);
atf=real(atf);
N1=round(Nt/3);

%Aceleração dada pela ode45
aode=(mu*(1-v.^2).*v -K.*(y+b.*y.^3))/M;

figure (4)
plot(t(N1:2*N1),aode(N1:2*N1),'r-',t(N1:2*N1),atf(N1:2*N1),'ob')
title('Aceleração')
xlabel('t(s)')
ylabel('a(t)')
legend('ode45','T.F.d')

%Diferença entre as acelerações obtidas 
Ndif=round(1024/3);
difa=abs(atf(Ndif:2*Ndif)-aode(Ndif:2*Ndif));
Dif_m=mean(difa)