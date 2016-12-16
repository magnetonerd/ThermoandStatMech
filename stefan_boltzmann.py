from numpy import *
from gaussquad import *

def f(z):
    return z**3/((1-z)**5*(exp(z/(1-z))-1))

N = 100
a = 0.0
b = 0.999
k = 1.3806e-23
T = 1.0
c = 2.9979e8
hbar = 1.0546e-34

const = k**4*T**4/(4*pi**2*c**2*hbar**3)

x,w = gaussxwab(N,a,b)

total = 0.0

for k in range(1,N):
    total += w[k]*f(x[k])

total = const*total
    
print total
