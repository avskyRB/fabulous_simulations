import numpy as np
import matplotlib.pyplot as plt

# This function improves the solution Ax=b with one Jacobi step, with the 
# result overwriting x. The matrix A is hardcoded and represents the Poisson equation.
    
def get_value(S, i, j):
    try:
        return S[i][j]
    except IndexError:
        return 0

def jacobi_step(x, b, N):
    xnew = np.zeros((N,N))
    #TODO: fill in your code
    for i in xrange(N):
    	for j in xrange(N):
    		xnew[i,j] =	 0.25 *(get_value(xnew,i-1,j)+get_value(xnew,i+1,j)+get_value(xnew,i,j-1)+get_value(xnew,i,j+1)-b[i,j])
    x[:,:] = xnew

# This function calculates the resdiuum vector res = b - Ax, for input vectors
# of length N. The output is stored in res.

'''
def gaus_siedel(x,b,N):
	xnew = np.zeros((N,N))
    for i in xrange(N):
    	for j in xrange(N):
    		xnew[i-1,j] = 0.25 *(get_value(xnew,i-2,j)+get_value(xnew,i,j)+get_value(xnew,i-1,j-1)+get_value(xnew,i-1,j+1)-b[i-1,j])
    		xnew[i,j-1] = 0.25 *(get_value(xnew,i-1,j-1)+get_value(xnew,i+1,j-1)+get_value(xnew,i,j-2)+get_value(xnew,i,j)-b[i,j-1])
    		xnew[i,j] =	 0.25 *(get_value(xnew,i-1,j)+get_value(xnew,i+1,j)+get_value(xnew,i,j-1)+get_value(xnew,i,j+1)-b[i,j])
    x[:,:] = xnew
'''


def calc_residuum(x, b, N, res):
    #TODO: fill in your code
    #res[i,j] = ...
    for i in xrange(N):
    	for j in xrange(N):
    		res[i,j] = b[i,j]+(get_value(x,i-1,j)+get_value(x,i+1,j)+get_value(x,i,j-1)+get_value(x,i,j+1)+4.0*get_value(x,i,j))/h*h
    
    
# This function calculates the norm of the vector of length N, 
# defined as the usual quadratic vector norm.   
def norm_of_residual(res, N):
    sum = np.sum(res**2)
    return np.sqrt(sum)
    
    
    
N = 256
steps = 2000
L = 1.0
h = L / N
eta = 0.1 * L
rho0 = 10.0

res = np.zeros((N,N))
  
#now set-up the density field
x = (np.arange(-N/2,N/2)+0.5)*h
mesh = np.meshgrid(x,x)
r2 = mesh[0] * mesh[0] + mesh[1] * mesh[1]
rho = rho0 * np.exp(-r2/(2*eta**2))
sum = np.sum(rho)

#initialize the starting values for phi[] and b[]
rho -= sum/(N**2)
b = 4*np.pi*h**2*rho
phi = np.zeros((N,N))

#open a file for outputting the residuum values, and then do 2000 Jacobi steps
residuals = np.zeros(steps)
f = file("res_jacobi.txt", "w")
f.write("%d\n"%steps)


for i in np.arange(steps):

	jacobi_step(phi, b, N)
	calc_residuum(phi, b, N, res)
	r = norm_of_residual(res, N)

	print("iter=%d:  residual=%g\n"%(i, r))
	f.write("%g\n"%r)
	residuals[i] = r
    
f.close()

plt.semilogy(residuals, label="jacobi")
plt.show()
