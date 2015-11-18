import numpy as np
import random

L = 1.0 # Length of the box
m = 1.0 # Mass of the particle
ngrid = 5 # Number of grids
h = L/ngrid
# Random position

pos = (0.5,0.5)
#pos = (random.uniform(0, L), random.uniform(0, L))

# DENSITY FIELD

#CIC assignment

rho = np.zeros((ngrid,ngrid)) # DENSITY
wp = np.zeros((ngrid,ngrid))

print pos

# FALTA TENER EN CUENTA QUE OCURRE EN LOS BORDES (A QUE CELDA VA LA CONTRIBUCION? IRIA A LA SIGUIENTE,
	#MISMO QUE DECIR QUE DE LA CAJA QUE HABRIA AL OTRO LADO ENTRARIA UNA HACIA ACA) BUSCAR MANERA INTELIGENTE DE PONER LOS IFS
for i in xrange(ngrid):
	for j in xrange(ngrid):
		if pos[0] - i*h < h and pos[0] - i*h > 0  and pos[1] - j*h < h and pos[1] - j*h > 0:
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

			



rho_kspace=np.fft.fft2(rho)

g = np.zeros((ngrid,ngrid))

def green(pos):
	for i in xrange(ngrid):
		for j in xrange(ngrid):
			g[i][j] = 1/(np.sqrt(((i-pos[0])*(i-pos[0])+(j-pos[1])*(j-pos[1]))))
	return g

green_real = green(pos)

green_kspace = np.fft.fft2(green_real)

potential_kspace = green_kspace*rho_kspace
potential = np.real(np.fft.ifft2(potential_kspace))

# OBTAINING THE FORCE

a = np.zeros((ngrid,ngrid,2))

for i in xrange(ngrid):
	for j in xrange(ngrid):
			print i,j
			if i != ngrid-1 and j != ngrid-1:
				a[i][j][0] = -(potential[i+1][j] - potential[i-1][j])/2*h
				a[i][j][1] = -(potential[i][j+1] - potential[i][j-1])/2*h
			if i== ngrid-1 and j != ngrid-1:
				a[i][j][0] = potential[i-1][j]/2*h
				a[i][j][1] = -(potential[i][j+1] - potential[i][j-1])/2*h

			if j==ngrid-1 and i != ngrid-1:
				a[i][j][0] = -(potential[i+1][j] - potential[i-1][j])/2*h
				a[i][j][1] = potential[i][j-1]/2*h

			if i== ngrid-1 and j == ngrid-1:
				a[i][j][0] = potential[i-1][j]/2*h
				a[i][j][1] = potential[i][j-1]/2*h



print potential
print a
