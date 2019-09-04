
% PROGRAMA DOS MATLAB MARIANA002
% ===========================================================
%% PROGRAMA DESARROLLADO POR RAYMUNDO JUÁREZ
% 25 DE AGOSTO DE 2017
% VERSIÓN 002 DEL FILTRO DE KALMAN EXTENDIDO
% PARA FUNCIÓN POLINÓMICA
% DE GRADO TRES

% LIMPIA VARIABLES, MEMORIA Y PANTALLA
% ADEMÁS CIERRA TODAS LAS VENTANAS ABIERTAS
clear all
clc 
close all

%% ASIGNACIÓN DE VALORES NUMÉRICOS COSNTANTES
T=0.01;

Q =0;
R = 1;
P0 = [];
x0 =0.99;

% FILTRO DE KALMAN EXTENDIDO
Qekf =0;
Rekf = 1;
P0ekf = [];
x0ekf =0.99;

P=1;
n = 36;

x = zeros(1,n);
z = zeros(n,1);

% RUIDO BLANCO
w = 0.5*randn(n,1);


%% SISTEMA
xs = @(xa) T+(1+T)*xa + (T/2)*xa^2+(T/6)*xa^3;
%xs = @(xa) p*m+(1+q-p)*xa - (q/m)*xa^2;

% CONDICIÓN INICIAL DEL SISTEMA
x(:,1) = x0;
z(1) = x0;

% ITERACIONES PARA CONSTRUIR EL SISTEMA
for k = 2:n
    x(:,k) = xs(x(:,k-1));
    z(k) = x(:,k) + w(k);
end

% GRAFICA DEL SISTEMA REAL
t = 0:T:(n-1)*T;
figure(1)
plot(t, x(1,:),'-r')
title('Real System')
%axis([0,0.35,0.5,4])
grid;
legend('True')


% INICIALIZACIÓN DE VARIABLES DEL FILTRO DE KALMAN EXTENDIDO
xpe = zeros(1,n);
xco = zeros(1,n);
Ppe = zeros(1,1,n);
Pco = zeros(1,1,n);


%CONDICION INICIAL PARA FILTRO KALMAN EXTENDIDO
xco0 = -10.09;
Pco0 = 1;
xco(:,1) = xco0;
Pco(:,:,1) = Pco0;

%FUNCIÓN DEL FILTRO DE KALMAN EXTENDIDO
xm = @(xd) T+(1+T)*xd + (T/2)*xd^2+(T/6)*xd^3;
Pm = @(xd,Pd) (1+T)^2*Pd + 2*(1+T)*T*xd*Pd + (1+2*T)*T*xd^2*Pd + ...
    T^2*xd^3*Pd+(T^2/4)*xd^4*Pd;

Kk = @(Pf) Pf*inv(Rekf+Pf); 
xM = @(xf,yf,Pf) xf+Kk(Pf)*(yf-xf);
PM = @(Pf) (1-Kk(Pf))^2 * Pf+Rekf*Kk(Pf)^2;

% ITERACIÓN PRINCIPAL DEL FILTRO
for k = 2:n
    % Prediccion Kalman Extendido
    xpe(:,k) = xm(xco(:,k-1));
    Ppe(:,:,k) = Pm(xco(:,k-1),Pco(:,:,k-1));
    
    % Correccion Kalman Extendido
    Kekf = Kk(Ppe(:,:,k));
    xco(:,k) = xM(xpe(:,k),z(k),Ppe(:,:,k));
    Pco(:,:,k) = PM(Ppe(:,:,k));
end

%GRÁFICA DE SISTEMA REAL Y ESTADO ESTIMADO POR FILTRO DE BASIN
t = 0:T:(n-1)*T;
figure(2)
plot(t, x(1,:),t, xco(1,:))
%axis([0,0.35,-15,5])
grid;
title('Estimation')
legend('True','Estimate EKF')

%GRÁFICA DE OBSERVACIONES
figure(3)
plot(t, z)
%axis([0,0.35,0.5,4])
grid;
title('Observations')
legend('Kalman Extended Observations')
