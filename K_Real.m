% PROGRAMA CINCO MATLAB K_REAL
% ======================================================
function dx=K_Real(T,x)
global A B C G V W L

dx = zeros(2,1); % un vector columna
X=[x(1);x(2)];
u=1;
dx=A*X + B*u +G*V*rand(1);
end