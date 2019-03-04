import numpy as np
import matplotlib.pyplot as plt
import time


#Matrix Lesson
A_c = np.array([[2,1,3,1],[2,6,8,3],[6,8,18,5]],dtype=float)

#Matrix ex 1
A_1 = np.array([[4,8,12,4],[3,8,13,5],[2,9,18,11]],dtype=float)

#Matrix ex 2
A_2 = np.array([[-27,9,-3,1,1],[-1,1,-1,1,-5],[8,4,2,1,3],[64,16,4,1,-5]],dtype=float)


###################### ######################
#                   Ex 1
###################### ######################

def PivotGauss(A):
	n = np.shape(A)[0]
	for k in range(0,n):
		piv = A[k,k]
		if(piv == 0):
			ind = np.argmax(A[:,k])
			temp = A[ind,:]
			A[ind,:] = A[k,:]
			A[k,:] = temp[:]
		if piv != 0:
			A[k,k:n+1] = A[k,k:n+1]/piv
			for i in range(k+1,n):
				A[i,k:n+1] = A[i,k:n+1] - A[i,k]*A[k,k:n+1]
				
	print(A)

	x = np.zeros(n)
	for i in range(0,n):
		ind = n-1-i
		x[ind] = A[ind,n] - sum(A[ind,ind+1:n]*x[ind+1:n])
		print ("x"+str(ind+1)+"=",x[ind])
	return x

print ("\nExercice 1 matrix solving :")
x = PivotGauss(A_1)


###################### ######################
#           Lesson matrix solution
###################### ######################

print ("\nLesson matrix solving :")
x = PivotGauss(A_c)

###################### ######################
#                   Ex 2
###################### ######################

print ("\nExercice 2 matrix solving :")
#Time estimation as in Exo 3
start = time.time()
x = PivotGauss(A_2)
end = time.time()

print("time :", end-start)
	
n = np.shape(A_2)[0]
def f(x,pol):
	somme = 0.
	for i in range(0,n):
		somme += pol[n-1-i]*x**i
	return somme


axis = np.linspace(-4.,5.,100)

data = np.zeros(100)
for i in range(0,100):
	data[i] = f(axis[i],x)


plt.plot(axis,data,linewidth=3)
plt.plot([-3,-1,2,4], [1,-5,3,-5], "o",markersize=8)
plt.show()


###################### ######################
#                   Ex 3
###################### ######################

#to be completed


