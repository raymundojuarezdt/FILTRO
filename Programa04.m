
% PROGRAMA CUATRO MATLAB KALMAN001
% ==========================================================
clear all
clc
global A B C G V W L Y n

n=1;
A=[-4 2;-2 -4]; B=[0;1]; C=[1,0]; G=[1;-1]; V=0.09; W=0.025;

[t,p]=ode45(@Ej_Kal,[0 10],[0.1 0 0.1]);

p=double(p);
S1=size(p);

for j=1:S1(1)
    
    P=[p(j,1) p(j,2);p(j,2) p(j,3)]
    L=P*C'*inv(W*rand(1))
    
    
    [T,x]=ode45(@K_Real,[0 10],[0.5 -0.5]);
    figure (1)
    plot(T,x(:,1),'--',T,x(:,2),'--')
    grid
    title('Real solution');
    xlabel('Time t');
    ylabel('Solution x');
    legend('x_1','x_2');
    
    S2=size(x);
    
    for i=1:S2(1)
        
        Y=C*[x(i,1);x(i,2)]+W*rand(1);
        
        
        [Te,xe]=ode45(@K_Estima,[0 10],[0.5 -0.5]);
        
    end
    
end
