import scipy.optimize

def f(x):
	return x**2 - 2.

def fprime(x):
	return 2.*x

print scipy.optimize.newton(f, 3., fprime)


def NewtonRaphson(f,fprime, x, eps=1e-8, max_iter=50):
	i = 0
	h = f(x)/fprime(x)
	while abs(h) >= eps and i <= max_iter:
		h = f(x)/fprime(x)
		x = x - h   
		#x(i+1) = x(i) - f(x)/f'(x)
		i += 1
	print x
	
NewtonRaphson(f,fprime,3.)
