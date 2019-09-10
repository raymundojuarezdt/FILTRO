using LinearAlgebra
using PyPlot
pygui(true)


#%%
F = [1 1; 0 1]
H = [1 0]
G = [0.5; 1]
g = 1

Q = [0 0;0 0]
R = 1

x0 = [100; 0]

n = 6
u = -g

x = zeros(2,n)
z = zeros(n,1)
w = 0.5*randn(n,1)

#%%
xs(xa) = F*xa + G*u'

#%%
x[:,1] = x0
for k = 2:n
    x[:,k] = xs(x[:,k-1])
    z[k] = H⋅x[:,k] + w[k] #OBSERVE: z recibe un escalar por lo que se usa dot() y no *
end

#%%
xp = zeros(2,n)
xc = zeros(2,n)

Pp = zeros(2,2,n)
Pc = zeros(2,2,n)

xc0 = [95; 1]
Pc0 = [10 0;0 1]

#%%
xc[:,1] = xc0
Pc[:,:,1] = Pc0

for k = 2:n
    # Prediccion
    xp[:,k]   = F*xc[:,k-1] + G*u
    Pp[:,:,k] = F*Pc[:,:,k-1]*F' + Q

    # Correccion
    S = H*Pp[:,:,k]⋅H' + R #OBSERVE: dot() en vez de *
    K = Pp[:,:,k]*H'*inv(S)
    xc[:,k] = xp[:,k] + K*(z[k] - H⋅xp[:,k]) #OBSERVE: dot() en vez de *
    Pc[:,:,k] = Pp[:,:,k] - K*S*K'
end

#%%
figure
t = 1:6
plot(t, x[1,:], t, z, t, xc[1,:])
axis([1,6,85,101])
grid(true)
legend(["True", "Measurement", "Estimate"])
