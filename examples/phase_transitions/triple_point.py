import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('triple_point.txt', delimiter='\s+')

def plot_graph(x, y):
	ax1 = df[:16].plot(kind='scatter', x=x, y=y, c='blue')
	df[16:25].plot(kind='scatter', x=x, y=y, c='red', ax=ax1)
	df[25:].plot(kind='scatter', x=x, y=y, c='black', ax=ax1)

plot_graph('T0', 'T')

plt.show()