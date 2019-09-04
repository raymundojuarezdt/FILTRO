% PROGRAMA TRES MATLAB MARIANA001
% % =======================================================
%% HOLA
% PROGRAMA DESARROLLADO POR RAYMUNDO JUÁREZ
% 27 DE JULIO DE 2017
% VERSIÓN 001 DEL FILTRO DE BASIN

clear all
clc 
close all

%% VALORES NUMÉRICOS
T=0.01;

% F = [1 1;0 1];
% H = [1 0];
% G = [0.5; 1];
% g = 1;
 Q =0;
 R = 1;
P0 = [];
x0 =0.99;

Qekf =0;
Rekf = 1;
P0ekf = [];
x0ekf =0.99;
P=1;
n = 36;
%u = -g;

x = zeros(1,n);
z = zeros(n,1);
w = 0.5*randn(n,1);


%SISTEMA ILTRO BASIN
xs = @(xa) T+(1+T)*xa + (T/2)*xa^2+(T/6)*xa^3;

%xs = @(xa) p*m+(1+q-p)*xa - (q/m)*xa^2;


%CONDICIÓN INICIAL
x(:,1) = x0;
z(1) = x0;


%ITERACIONES PARA CONSTRUIR EL SISTEMA
for k = 2:n
    x(:,k) = xs(x(:,k-1));
    z(k) = x(:,k) + w(k);
    
    
end

%GRAFICAR EL SISTEMA REAL
t = 0:T:(n-1)*T;
figure(1)
plot(t, x(1,:),'-r')
title('Real System')
axis([0,0.35,0.5,4])
grid;
legend('True Basin')


%VARIABLES DEL FILTRO DE BASIN
xp = zeros(1,n);
xc = zeros(1,n);
Pp = zeros(1,1,n);
Pc = zeros(1,1,n);

xpe = zeros(1,n);
xco = zeros(1,n);
Ppe = zeros(1,1,n);
Pco = zeros(1,1,n);

%CONDICION INICIAL PARA FILTRO BASIN
xc0 = -10.09;
Pc0 = 1;
xc(:,1) = xc0;
Pc(:,:,1) = Pc0;

xco0 = -10.09;
Pco0 = 1;
xco(:,1) = xco0;
Pco(:,:,1) = Pco0;

%FUNCIÓN DEL FILTRO DE BASIN
xe = @(xb,Pb) T+(1+T)*xb + (T/2)*xb^2+(T/6)*xb^3+(T/2)*Pb+(T/2)*xb*Pb;
Pe = @(xb,Pb) (1+T)^2*Pb + 2*(1+T)*T*xb*Pb+(1+3*T/2)*T*Pb^2 + (1+2*T)*T*xb^2*Pb + ...
    (5*T^2/12)*Pb^3 + T^2*xb^2*Pb^2 + (2*T^2/9)*xb^4*Pb + (2/3)*(1+T)*T*xb^4 + ...
    2*T^2*xb*Pb^2 + T^2*xb^3*Pb;

K = @(Pc) Pc*inv(R+Pc); 
xo = @(xc,y,Pc) xc+K(Pc)*(y-xc);
Po = @(Pc) (1-K(Pc))^2 * Pc+R*K(Pc)^2;

%FUNCIÓN DEL FILTRO DE EKF
xm = @(xd) T+(1+T)*xd + (T/2)*xd^2+(T/6)*xd^3;
Pm = @(xd,Pd) (1+T)^2*Pd + 2*(1+T)*T*xd*Pd + (1+2*T)*T*xd^2*Pd + ...
    T^2*xd^3*Pd+(T^2/4)*xd^4*Pd;

Kk = @(Pf) Pf*inv(Rekf+Pf); 
xM = @(xf,yf,Pf) xf+K(Pf)*(yf-xf);
PM = @(Pf) (1-K(Pf))^2 * Pf+Rekf*K(Pf)^2;



%ITERACIÓN PRINCIPAL DEL FILTRO
for k = 2:n
    % Prediccion
    xp(:,k) = xe(xc(:,k-1),Pc(:,:,k-1));
    Pp(:,:,k) = Pe(xc(:,k-1),Pc(:,:,k-1));
    
    xpe(:,k) = xm(xco(:,k-1));
    Ppe(:,:,k) = Pm(xco(:,k-1),Pco(:,:,k-1));
    
    % Correccion
    
    K1 = K(Pp(:,:,k));
    xc(:,k) = xo(xp(:,k),z(k),Pp(:,:,k));
    Pc(:,:,k) = Po(Pp(:,:,k));
    
    Kekf = Kk(Ppe(:,:,k));
    xco(:,k) = xM(xpe(:,k),z(k),Ppe(:,:,k));
    Pco(:,:,k) = PM(Ppe(:,:,k));
end
%grafica finales
t = 0:T:(n-1)*T;
figure(2)
plot(t, x(1,:),t, xc(1,:),t, xco(1,:))
axis([0,0.35,-15,5])
grid;
title('Estimation')
legend('True', 'Estimate','Estimate EKF')

figure(3)
plot(t, z)
axis([0,0.35,0.5,4])
grid;
title('Observations')
legend('Basin Obs')
