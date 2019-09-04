% PROGRAMA SEIS MATLAB K_ESTIMA
% ==========================================================
function dxe=K_Estima(Te,xe)
global A B C G V W L Y

dxe = zeros(2,1); % un vector columna
X=[xe(1);xe(2)];
u=1;
dxe=A*X + B*u +L*(Y-C*X);
end