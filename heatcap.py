"""
This program computes the volume heat capacity
for solids based on the Debye Theory. It then
Plots this against different temperatures.
"""
from pylab import plot,show
from numpy import *

#The next three functions are gaussian quadrature. Not important part of code
#For full explanation go to github.com/magnetonerd/mathlib/gaussquad.py

#######################DON'T FOCUS ON THIS PART OF THE CODE####################
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
######################IMPORTANT PART OF CODE###########################
"""
The following solve the integral
                                ___ThetaD/T
                                |               x^4 * exp(x)
Cv = 9*Vol*rho*kB*(T/ThetaD)^3  |              ______________ dx
                                |               (exp(x) - 1)^2
                              ___ 0.0
"""
Vol    = 0.001      #Volume of solid             m^3
rho    = 6.022e28   #Number density of particles 1/m^3
thetaD = 428.0      #Debye temperature           K
kB     = 1.3806e-23 #Boltzman's Constant         J/K
Temp   = [0.0]*495  #Temperature of solid        K
Cv     = [0.0]*495  #Volume Heat Capacity        J/(kg K)

for T in range(5,500):
    C = 9*Vol*rho*kB*T**3/thetaD**3
    N = 100 #Number of sample point for approximating the integral
    a = 0.0 #Initial value of integral
    b = thetaD/T #final value of integral
    x,w = gaussxwab(N,a,b) #Getting the sample points and weights
    total = 0.0 #initializing the value of the integral
    #Performing the integral
    for k in range(1,N):
        total += w[k]*f(x[k])
    total = total*C#multiplying by the constants
    Cv[T-5],Temp[T-5] = total,T #setting up arrays to facilitate plotting
    
plot(Temp,Cv) #generating the plot
show() #Making sure people can see the plot.
