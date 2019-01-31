import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.integrate as integ

x = np.linspace(0,5,100)
y = sp.sin(x) + sp.cos(2*x)

print integ.simps(y, x)

def f_int(x):
	return -sp.cos(x) + sp.sin(2*x)/2
	
print (f_int(5) - f_int(0))	

plt.plot(x,y)
plt.show()
