from adams import Integrator
from general import Satellite
from math import sqrt, pi

# Crating integrator object
I = Integrator(3)

# Calculations of orbit parameters
Mu = 3.986004415e14
e = 0.0001
x = (6380 + 700)*1e3
a = x / (1 - e)
V = sqrt(Mu / a * (1 + e) / (1 - e))
T = 2 * pi * sqrt(a**3 / Mu) * 10

s = Satellite(280, x, 0, 0, 0, 0, V)

# Starting simulation by using chosen integrator
I.integrate(s, T)
