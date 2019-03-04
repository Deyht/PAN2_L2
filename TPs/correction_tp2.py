import numpy as np
import matplotlib.pyplot as plt

#Problem parameters initialization
m = 5.0
g = 9.81
v0 = 30
theta = np.pi/4.
ti = 0
tf = 5
dt = 0.1

t = np.arange(0,tf,dt) 
#Warning : wille give N points -> np.arange(0,tf+dt,dt) for N+1 points 
N = int(tf/dt) # The number of steps MUST be integer



######################## ########################
#               solution analytique
######################## ########################
a_r = np.zeros((2,N),dtype = float)
a_v = np.zeros((2,N),dtype = float)
a_r[:,0] = 0.0
vx_0 = v0*np.cos(theta)
vy_0 = v0*np.sin(theta)

for i in range(0,N):
	a_r[0,i] = vx_0*(i*dt)
	a_r[1,i] = -0.5*g*(i*dt)**2 + vy_0*(i*dt)
	a_v[0,i] = vx_0 
	a_v[1,i] = -g*(i*dt) + vy_0




######################## ########################
#            Fonction second membre
######################## ########################

def f(r, v):
	vect = np.zeros(4)
	vect[0] = v[0]
	vect[1] = v[1]
	vect[2] = 0.0
	vect[3] = -g
	return vect
	
	
	

######################## ########################
#                Methode d'euler
######################## ########################
r = np.zeros((2,N),dtype = float)
v = np.zeros((2,N),dtype = float)

r[:,0] = 0.0
v[0,0] = vx_0
v[1,0] = vy_0

for i in range(0,N-1):
	fc = f(r[:,i],v[:,i])
	r[:,i+1] = r[:,i] + fc[0:2]*dt
	v[:,i+1] = v[:,i] + fc[2:4]*dt
	
	
	
	
	
######################## ########################
#                  Methode RK4
######################## ########################
r2 = np.zeros((2,N),dtype = float)
v2 = np.zeros((2,N),dtype = float)

r2[:,0] = 0.0
v2[0,0] = vx_0
v2[1,0] = vy_0

for i in range(0,N-1):
	k1 = f(r2[:,i],v2[:,i])
	k2 = f(r2[:,i]+ k1[0:2]*dt*0.5,v2[:,i]+ k1[2:4]*dt*0.5)
	k3 = f(r2[:,i]+ k2[0:2]*dt*0.5,v2[:,i]+ k2[2:4]*dt*0.5)
	k4 = f(r2[:,i]+ k3[0:2]*dt,v2[:,i]+ k3[2:4]*dt)
	k = (k1 + 2.*k2 + 2.*k3 + k4)/6.
	r2[:,i+1] = r2[:,i] + dt*k[0:2]
	v2[:,i+1] = v2[:,i] + dt*k[2:4]
	
	
	
	
	
######################## ########################
#             Methode Velocity-Verlet
######################## ########################

r3 = np.zeros((2,N),dtype = float)
v3 = np.zeros((2,N),dtype = float)

r3[:,0] = 0.0
v3[0,0] = vx_0
v3[1,0] = vy_0

a = np.array([0,-g])

for i in range(0,N-1):
	r3[:,i+1] = r3[:,i] + v3[:,i]*dt + 0.5*(dt**2)*a[:]
	#In our case acceleraction is constant over time, so da/dt = 0
	v3[:,i+1] = v3[:,i] + a[:]*dt
	





######################## ########################
#				 Ploting solutions
######################## ########################

###### CHANGE COMMENTED LINES TO SEE CHANGES BETWEEN METHODS ######
plt.plot(a_r[0,:], a_r[1,:], label="analytique")
plt.plot(r[0,:],r[1,:], label="Euler")
plt.plot(r2[0,:],r2[1,:],label="RK4")
plt.plot(r3[0,:],r3[1,:],label="Velocity-Verlet")
plt.legend()
plt.savefig("pos.png")

plt.clf()

###### CHANGE COMMENTED LINES TO SEE CHANGES BETWEEN METHODS ######
## POSITION ERRORS ##
#plt.plot(t[:], a_r[0,:]-r[0,:],label="E_x EULER")
#plt.plot(t[:], a_r[1,:]-r[1,:],label="E_y EULER")
#plt.plot(t[:], a_r[0,:]-r2[0,:],label="E_x RK4")
plt.plot(t[:], a_r[1,:]-r2[1,:],label="E_y RK4")
#plt.plot(t[:], a_r[0,:]-r3[0,:],label="E_x V-V")
plt.plot(t[:], a_r[1,:]-r3[1,:],label="E_y V-V")

## VELOCITY ERRORS ##
#plt.plot(t[:], a_v[0,:]-v[0,:],label="E_vx EULER")
#plt.plot(t[:], a_v[1,:]-v[1,:],label="E_vy EULER")
#plt.plot(t[:], a_v[0,:]-v2[0,:],label="E_vx RK4")
#plt.plot(t[:], a_v[1,:]-v2[1,:],label="E_vy RK4")
#plt.plot(t[:], a_v[0,:]-v3[0,:],label="E_vx V-V")
#plt.plot(t[:], a_v[1,:]-v3[1,:],label="E_vy V-V")
plt.legend()
plt.savefig("error.png")

plt.clf()

######################## ######################## ########################




######################## ########################
#				       EX2
######################## ########################

from scipy.integrate import odeint

pos = np.zeros((2,N),dtype = float)

dt = 0.1
tf = 30.
N = int(tf/dt)

t = np.arange(0,tf,dt)
#tab for euler
pos = np.zeros((2,N),dtype = float)
pos[0,0] = 1.
pos[1,0] = 0.

#tab for RK4
pos2 = np.zeros((2,N),dtype = float)
pos2[0,0] = 1.
pos2[1,0] = 0.

#Define the seconf member function f based on the EDO system definition
def derivee(vec,t):
	return np.array([-0.25*vec[0] + vec[1], -vec[0] - 0.25*vec[1]])

#Solving using Euler
for i in range(0,N-1):
	f = derivee(pos[:,i],i*dt)
	pos[0,i+1] = pos[0,i] + f[0]*dt
	pos[1,i+1] = pos[1,i] + f[1]*dt


#Solving using RK4
for i in range(0,N-1):
	k1 = derivee(pos2[:,i],i*dt)
	k2 = derivee(pos2[:,i] + k1[:]*dt*0.5,i*dt)
	k3 = derivee(pos2[:,i] + k2[:]*dt*0.5,i*dt)
	k4 = derivee(pos2[:,i] + k3[:]*dt,i*dt)
	k = (k1 + 2.*k2 + 2.*k3 + k4)/6.
	pos2[:,i+1] = pos2[:,i] + dt*k[:]


#scipy.integrate.odeint(func, y0, t, ...)
# func = fonction to compute the derivative
# y0 = inital conditions
# t = time points array
#return the corresponding tab, based on t, of the same size
ode_pos = odeint(derivee,[1,0],t)


plt.plot(pos[0,:],pos[1,:], label="Euler")
#plt.plot(pos2[0,:],pos2[1,:], label="RK4")
plt.plot(ode_pos[:,0], ode_pos[:,1], label="odeint")
plt.legend()
plt.savefig("ex2_pos.png")








