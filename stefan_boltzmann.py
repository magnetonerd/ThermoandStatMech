from numpy import * #exp
from gaussquad import * #gaussxwab

"""
This program calculates the stefan-boltzmann constant.
It does this by using the gaussian quadrature method.
This is because the function is a relatively smooth
function.
"""

"""
The original function is         
                                  x^3
                              __________ dx
                              exp(x) - 1
However this is integrated over an infinite range.
So, there was a change of variables:
                                 z
                         x =  _______
                               1 - z
Which maps from [0,inf) to [0,1]
"""

def f(z):
    return z**3/((1-z)**5*(exp(z/(1-z))-1))

"""
The Energy function is
                       ___ 1
         k^4 T^4      |
W = _________________ |        z^3             1
    4 pi^2 c^2 hbar^3 |     _________  _________________  dz = sigma T^4
                      |     (1 - z)^5  (exp(z/(1-z)) - 1
                    ___ 0
"""

N = 100           # Number of steps in integration
a = 0.0           # lower bound
b = 0.999         # just to stop an overflow warning
k = 1.3806e-23    # Boltzmann's Constant J/K
c = 2.9979e8      # Speed of Light m/s
hbar = 1.0546e-34 # Reduced Planck's Constant Js

const = k**4/(4*pi**2*c**2*hbar**3) # Combining all the constants into a single constant

x,w = gaussxwab(N,a,b) #Getting the sample points and weights for integration

total = 0.0 #initializing the integral

#Calculating the integral
for k in range(1,N):
    total += w[k]*f(x[k])

total = const*total #multiplying the integral by all the constants
    
print total #displaying the result

"""
The result works out to be 5.6692079303e-8.
The real value is          5.670367(13)e-8.
The agreement can be improved by playing
with precision of the constants or with the
number of steps in the integration.
"""
