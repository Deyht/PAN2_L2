import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint


G = 6.67408e-11
M = 1.98855e30
UA = 1.496e11
a = UA
e = 0.0

nb_j_a = 365.24218979
dt = 0.02*24*3600
dt_o = dt
periode = 2*np.pi*np.sqrt(a**3/(G*M))#/(24.0*3600.)
tf = 2*periode
N = int(tf/dt)
print N

fact = G*M/UA**3
rot = 2*np.pi/periode
print rot
print fact
print periode

#solution analytique
t = np.arange(0,tf,dt)
pos_a = np.zeros((2,N+1), dtype=float)

pos_a[0,:] = a*np.cos(rot*t[:])
pos_a[1,:] = a*np.sin(rot*t[:])


pos_1 = np.zeros((2,N+1), dtype=float)
vit_1 = np.zeros((2,N+1), dtype=float)
pos_1[:,0] = (a,0)
vit_1[:,0] = (0,2*np.pi*a/periode)

#Methode d'euler
for i in range (0, N):
	r = np.sqrt(pos_1[0,i]**2 + pos_1[1,i]**2)
	pos_1[:,i+1] = pos_1[:,i] + vit_1[:,i]*dt
	vit_1[:,i+1] = vit_1[:,i] - (G*M/r**3)*pos_1[:,i]*dt


pos_2 = np.zeros((2,N+1), dtype=float)
vit_2 = np.zeros((2,N+1), dtype=float)
pos_2[:,0] = (a,0)
vit_2[:,0] = (0,2*np.pi*a/periode)

#methode Velocity-Verlet
for i in range (0,N):
	r = np.sqrt(pos_2[0,i]**2 + pos_2[1,i]**2)
	pos_2[:,i+1] = pos_2[:,i] + vit_2[:,i]*dt - 0.5*(dt**2)*(G*M/r**3)*pos_2[:,i]
	#In our case acceleraction is constant over time, so da/dt = 0
	vit_2[:,i+1] = vit_2[:,i] - (G*M/r**3)*pos_2[:,i]*dt

#methode odeint
def derivee(vec,t):
	r = np.sqrt(vec[0]**2 + vec[1]**2)
	return np.array([vec[2], vec[3], -(G*M/r**3)*vec[0] , -(G*M/r**3)*vec[1]])

odeint_pos = odeint(derivee,[a,0.0,0.0,2*np.pi*a/periode],t)


plt.plot(pos_2[0,:], pos_2[1,:], label="V-V")
plt.plot(pos_1[0,:], pos_1[1,:], label="Euler")
plt.plot(pos_a[0,:], pos_a[1,:], label="Analytic")
plt.plot(odeint_pos[:,0], odeint_pos[:,1], label="Odeint") 
plt.legend()
plt.title("Position plot using SI with dt ="+str(dt_o), y = 1.08)
plt.savefig("pos_si"+str(dt_o)+".pdf")
plt.clf()

r_a = np.sqrt(pos_a[0,:]**2 + pos_a[1,:]**2)
r_1 = np.sqrt(pos_1[0,:]**2 + pos_1[1,:]**2)
r_2 = np.sqrt(pos_2[0,:]**2 + pos_2[1,:]**2)
r_o = np.sqrt(odeint_pos[:,0]**2 + odeint_pos[:,1]**2)

plt.plot(t, (r_o - r_a)/r_a, label="Error Odeint")
plt.plot(t, (r_1 - r_a)/r_a, label="Error Euler")
plt.plot(t, (r_2 - r_a)/r_a, label="Error V-V")
plt.legend()
plt.title("Error plot using SI with dt ="+str(dt_o), y = 1.08)
plt.savefig("err_si"+str(dt_o)+".pdf")
plt.clf()




###################### ######################
#             With proper units
###################### ######################


dt = dt/periode
tf = 2
pos_p = np.zeros((2,N+1),dtype=float)
vit_p = np.zeros((2,N+1),dtype=float)
pos_p[:,0] = (1.0,0)
vit_p[:,0] = (0,2*np.pi)

t = np.arange(0,tf,dt)

pos_pa = np.zeros((2,N+1), dtype=float)

pos_pa[0,:] = np.cos(rot*t[:])
pos_pa[1,:] = np.sin(rot*t[:])

#methode Velocity-Verlet
for i in range (0,N):
	r = np.sqrt(pos_p[0,i]**2 + pos_p[1,i]**2)
	pos_p[:,i+1] = pos_p[:,i] + vit_p[:,i]*dt - 0.5*(dt**2)*(2*np.pi)**2/(r**3)*pos_p[:,i]
	#In our case acceleraction is constant over time, so da/dt = 0
	vit_p[:,i+1] = vit_p[:,i] - (2*np.pi)**2/(r**3)*pos_p[:,i]*dt


pos_p2 = np.zeros((2,N+1),dtype=float)
vit_p2 = np.zeros((2,N+1),dtype=float)
pos_p2[:,0] = (1.0,0)
vit_p2[:,0] = (0,2*np.pi)

#Methode d'Euler
for i in range (0, N):
	r = np.sqrt(pos_p2[0,i]**2 + pos_p2[1,i]**2)
	pos_p2[:,i+1] = pos_p2[:,i] + vit_p2[:,i]*dt
	vit_p2[:,i+1] = vit_p2[:,i] - (2*np.pi)**2/(r**3)*pos_p2[:,i]*dt

#methode odeint
def derivee2(vec,t):
	r = np.sqrt(vec[0]**2 + vec[1]**2)
	return np.array([vec[2], vec[3], -(2*np.pi)**2/(r**3)*vec[0] , -(2*np.pi)**2/(r**3)*vec[1]])


odeint_pos2 = odeint(derivee2,[1,0,0,2*np.pi],t/periode)

plt.plot(odeint_pos2[:,0], odeint_pos2[:,1], label="Odeint")
plt.plot(pos_p2[0,:], pos_p2[1,:], label="Euler")
plt.plot(pos_p[0,:], pos_p[1,:], label = "V-V")
plt.plot(pos_pa[0,:], pos_pa[1,:], label="Analytic")
plt.legend()
plt.title("Position plot using NAT with dt ="+str(dt_o), y = 1.08)
plt.savefig("pos_nat"+str(dt_o)+".pdf")
plt.clf()

r_pa = np.sqrt(pos_pa[0,:]**2 + pos_pa[1,:]**2)
r_p = np.sqrt(pos_p[0,:]**2 + pos_p[1,:]**2)
r_p2 = np.sqrt(pos_p2[0,:]**2 + pos_p2[1,:]**2)
r_o2 = np.sqrt(odeint_pos2[:,0]**2 + odeint_pos2[:,1]**2)

plt.plot(t, (r_p - r_pa)/r_pa, label="Error Euler")
plt.plot(t, (r_p2 - r_pa)/r_pa, label="Error V-V")
plt.plot(t, (r_o2 - r_pa)/r_pa, label="Error Odeint")
print(np.max(np.abs((r_p - r_pa)/r_pa)),np.max(np.abs((r_p2 - r_pa)/r_pa)),np.max(np.abs((r_o2 - r_pa)/r_pa)))
plt.legend()
plt.title("Error plot using NAT with dt ="+str(dt_o), y = 1.08)
plt.savefig("err_nat"+str(dt_o)+".pdf")
plt.clf()




###################### ######################
#                   Ex 2
###################### ######################

nb_j_a = 365.24218979
dt = 0.001*24*3600 #jours terrestres
dt_o = dt
a = 0.387
e = 0.206
alpha = 0.001
periode_m = 2*np.pi*np.sqrt((a*UA)**3/(G*M))
tf = 5*periode_m/periode


#Periode terrestre 
dt = dt/periode
t = np.arange(0,tf,dt)
N = int(tf/dt)

pos_m = np.zeros((2,N+1),dtype=float)
vit_m = np.zeros((2,N+1),dtype=float)
pos_m[:,0] = ((1+e)*a,0)
vit_m[:,0] = (0,2*np.pi*np.sqrt((1-e)/(a*(1+e))))


#methode Velocity-Verlet
for i in range (0,N):
	r = np.sqrt(pos_m[0,i]**2 + pos_m[1,i]**2)
	pos_m[:,i+1] = pos_m[:,i] + vit_m[:,i]*dt - 0.5*(dt**2)*(2*np.pi)**2/(r**3)*pos_m[:,i]
	#In our case acceleraction is constant over time, so da/dt = 0
	vit_m[:,i+1] = vit_m[:,i] - (2*np.pi)**2/(r**3)*pos_m[:,i]*dt


plt.plot(pos_m[0,:], pos_m[1,:])
plt.savefig("pos_mercure_VV_"+str(dt_o)+".pdf")
plt.clf()



###################### ######################
#                   Ex 3
###################### ######################

pos_m2 = np.zeros((2,N+1),dtype=float)
vit_m2 = np.zeros((2,N+1),dtype=float)
pos_m2[:,0] = ((1+e)*a,0)
vit_m2[:,0] = (0,2*np.pi*np.sqrt((1-e)/(a*(1+e))))


#methode Velocity-Verlet
for i in range (0,N):
	r = np.sqrt(pos_m2[0,i]**2 + pos_m2[1,i]**2)
	pos_m2[:,i+1] = pos_m2[:,i] + vit_m2[:,i]*dt - 0.5*(dt**2)*(2*np.pi)**2/(r**3)*(1+alpha/r**2)*pos_m2[:,i]
	#In our case acceleraction is constant over time, so da/dt = 0
	vit_m2[:,i+1] = vit_m2[:,i] - (2*np.pi)**2/(r**3)*(1+alpha/r**2)*pos_m2[:,i]*dt


plt.plot(pos_m2[0,:], pos_m2[1,:])
plt.savefig("pos_mercure_precess_VV_"+str(dt_o)+".pdf")
plt.clf()


r_m = np.sqrt(pos_m2[0,:]**2 + pos_m2[1,:]**2)
dr_m = (pos_m2[0,:]*vit_m2[0,:] + pos_m2[1,:]*vit_m2[1,:])/r_m[:]
theta_m = np.arctan2(pos_m2[1,:], pos_m2[0,:])

print np.shape(r_m), np.shape(dr_m), np.shape(theta_m)

plt.plot(t, theta_m, label="Angle")
plt.plot(t, r_m, label="r")
plt.plot(t, dr_m, label="dr")
plt.legend()
plt.savefig("orb_param_mercure"+str(dt_o)+".pdf")
plt.clf()


theta_aph = []
time = []
for i in range(0,N):
	if((dr_m[i+1]/dr_m[i] < 0.0) & (theta_m[i] >= 0.)):
		theta_aph.append(theta_m[i])
		time.append(i*dt)

print theta_aph

#for true aphelie value take only positiv absolute value
plt.plot(time, theta_aph)
plt.savefig("precess_mercure_"+str(dt_o)+".pdf")

deriv = ((theta_aph[0] - theta_aph[-1])/(time[0] - time[-1]))

deriv *= ((180.0/(np.pi))*3600.)*(100.0)*(1.1e-8/alpha)

#value of perihelion precess un arc sec per century (observed = 43''/century)
print deriv













