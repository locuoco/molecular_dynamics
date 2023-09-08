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

print('calculated T_c =', df['T'][20])
print('calculated P_c =', df['P'][20])

prop_cst = eps / boltzmann
print('Schultz T_c =', 1.321 * prop_cst, '+/-', 0.007 * prop_cst)

prop_cst = eps / sigma**3 / 101325
print('Schultz P_c =', 0.129 * prop_cst, '+/-', 0.005 * prop_cst)

def plot_graph(x, y):
	ax1 = df[:20].plot(kind='scatter', x=x, y=y, c='black')
	df[20:].plot(kind='scatter', x=x, y=y, c='red', ax=ax1)

plot_graph('T0', 'T')
plot_graph('T0', 'P')
plot_graph('T', 'P')

plt.show()
































