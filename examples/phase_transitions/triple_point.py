import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

eps = 1.54e-21
sigma = 3.623e-10
boltzmann = 1.380649e-23

df = pd.read_csv('triple_point.txt', delimiter='\s+')

def plot_graph(x, y):
	ax1 = df[:16].plot(kind='scatter', x=x, y=y, c='blue')
	df[16:25].plot(kind='scatter', x=x, y=y, c='red', ax=ax1)
	df[25:].plot(kind='scatter', x=x, y=y, c='black', ax=ax1)

plot_graph('T0', 'T')
#plot_graph('P', 'T')

crit_temps = df['T'][16:25]

prop_cst = eps / boltzmann

print('calculated T3 =', np.mean(crit_temps), '+/-', np.std(crit_temps) / np.sqrt(len(crit_temps)-1))
print('Schultz T3 =', 0.690 * prop_cst, '+/-', 0.005 * prop_cst)

plt.show()