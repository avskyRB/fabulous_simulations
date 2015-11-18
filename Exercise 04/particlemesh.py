import numpy as np
import random

L = 1.0 # Length of the box
m = 1.0 # Mass of the particle
ngrid = 5 # Number of grids
h = L/ngrid
# Random position

#pos = (random.uniform(0, L), random.uniform(0, L))
pos=(0.5,0.5)

# DENSITY FIELD

#CIC assignment

rho = np.zeros((ngrid,ngrid)) # DENSITY
wp = np.zeros((ngrid,ngrid))

print pos

# FALTA TENER EN CUENTA QUE OCURRE EN LOS BORDES (A QUE CELDA VA LA CONTRIBUCION? IRIA A LA SIGUIENTE,
	#MISMO QUE DECIR QUE DE LA CAJA QUE HABRIA AL OTRO LADO ENTRARIA UNA HACIA ACA) BUSCAR MANERA INTELIGENTE DE PONER LOS IFS
	print NO MORE SPANISH COMMENT
for i in xrange(ngrid):
	for j in xrange(ngrid):
		if pos[0] - i*h < h and pos[0] - i*h > 0  and pos[1] - j*h < h and pos[1] - j*h > 0:
			print i,j
			pf = pos[0]/h - 0.5
			p = i
			pa = pf - p
			qf = pos[1]/h -0.5
			q = j
			qa = qf - q
			wp[i][j] = (1-pa)*(1-qa)
			rho[i][j] =(1/h*h)*(m*wp[i][j])
			if i != ngrid-1:
				wp[i+1][j] = pa * (1-qa)
				rho[i+1][j] =(1/h*h)*(m*wp[i+1][j])
			if i == ngrid-1:
				wp[0][j] = pa * (1-qa)
				rho[0][j] =(1/h*h)*(m*wp[0][j])
			if j != ngrid-1:
				wp[i][j+1] = (1-pa)*qa
				rho[i][j+1] =(1/h*h)*(m*wp[i][j+1])
			if j == ngrid-1:
				wp[i][0] = (1-pa)*qa
				rho[i][0] =(1/h*h)*(m*wp[i][0])
			if i != ngrid-1 and j!= ngrid-1:
				wp [i+1][j+1] = pa*qa
				rho[i+1][j+1] =(1/h*h)*(m*wp[i+1][j+1])
			if i == ngrid-1 and j == ngrid-1:
				wp[0][0] = pa*qa 
				rho[0][0] =(1/h*h)*(m*wp[0][0])

			


#print wp
#print rho

rho_kspace=np.fft.fft2(rho)
green = np.zeros((ngrid,ngrid))

for i in xrange(ngrid):
    for j in xrange(ngrid):
        green[i][j]=1/np.sqrt(np.power(pos[0]-i,2)+np.power(pos[1]-j,2)) #green function in matrix form
        
print green
        
green_kspace=np.fft.fft2(green) #fourier transform of the green function
phi_kspace=green_kspace*rho_kspace #convolution product
phi_real=np.real(np.fft.ifft2(phi_kspace)) #potential

print "          "

print phi_real
