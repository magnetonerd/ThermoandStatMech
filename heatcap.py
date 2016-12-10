from pylab import plot,show
from numpy import *

EPS = 1e-15

def gaussxw(N):
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a)
    
    delta = 1.0
    while delta>EPS:
        p0 = [1.0]*N
        p1 = x
        for k in range(1,N):
            p0,p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dP = (N+1)*(x*p1-p0)/(x**2-1)
        dx = p1/dP
        x -= dx
        delta = max(abs(dx))

    w = 2*(N+1)**2/(N**2*(1-x**2)*dP**2)
    return x,w

def gaussxwab(N,a,b):
    x,w=gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

def f(x):
    return x**4*exp(x)/(exp(x)-1)**2

Vol    = 0.001      #m^3
rho    = 6.022e28   #1/m^3
thetaD = 428.0      #K
kB     = 1.3806e-23 #J/K
Temp   = [0.0]*495  #K
Cv     = [0.0]*495  #J/(kg K)

for T in range(5,500):
    C = 9*Vol*rho*kB*T**3/thetaD**3
    N = 100
    a = 0.0
    b = thetaD/T
    x,w = gaussxwab(N,a,b)
    total = 0.0
    for k in range(1,N):
        total += w[k]*f(x[k])
    total = total*C
    Cv[T-5],Temp[T-5] = total,T
    
plot(Temp,Cv)
show()
