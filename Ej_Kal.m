% PROGRAMA OCHO MATLAB EJ_KAL
% ===========================================================

function dp=Ej_Kal(t,p)
global A B C G V W L

dp = zeros(3,1); % un vector columna
P=[p(1),p(2);p(2),p(3)];
DP=A*P + P*A' - P*C'*inv(W*rand(1))*C*P + G*V*rand(1)*G';
dp(1)=DP(1,1);
dp(2)=DP(1,2);
dp(3)=DP(2,2);
end
