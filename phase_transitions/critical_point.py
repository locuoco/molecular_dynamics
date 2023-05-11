import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

eps = 1.54e-21
sigma = 3.623e-10
boltzmann = 1.380649e-23

def lofti_vaporliquid(temperature):
	"""
	Computes the vapor-liquid coexistence pressure for a Lennard-Jones substance
	as a function of temperature using Lofti approximation (1992)

	  temperature: temperature in kelvin (K)
	  return: pressure in atmospheres (atm)
	"""
	reduced_temp = boltzmann * temperature / eps
	reduced_pressure = np.exp(1.2629*reduced_temp - 4.8095/reduced_temp - 0.15115/reduced_temp**4)
	pressure = reduced_pressure * eps / sigma**3
	return pressure / 101325

df = pd.read_csv('critical_point.txt', delimiter='\s+')

df['P'] /= 1.458397257425908222e-5

print('T =', df['T'][20])
print('P =', df['P'][20])

df.plot(kind='scatter', x='T0', y='T', c='blue')
df.plot(kind='scatter', x='T0', y='P', c='blue')

#p = np.arange(120, 150.8, 1)

fig, ax1 = plt.subplots(1)
#ax1.plot(p, lofti_vaporliquid(p), c='black')
df.plot(kind='scatter', x='T', y='P', c='blue', ax=ax1)

plt.show()
































