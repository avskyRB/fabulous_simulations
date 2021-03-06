import numpy as np
import random

L = 1.0 # Length of the box
m = 1.0 # Mass of the particle
ngrid = 5 # Number of grids
h = L/ngrid
# Random position

#pos = (0.5,0.5)
pos = (random.uniform(0, L), random.uniform(0, L))

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

a = np.zeros((ngrid,ngrid,2)) # Aceleration field (matrix NxNx2)

for i in xrange(ngrid):
	for j in xrange(ngrid):
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



#CIC interpolation of the force
F = np.zeros((ngrid,ngrid,2))

for i in xrange(ngrid):
	for j in xrange(ngrid):
		F[i][j][0] = m*a[i][j][0]*wp[i][j]
		F[i][j][1] = m*a[i][j][1]*wp[i][j]


# distances distributed around:
r_min = 0.3*L/ngrid
r_max=L/2
p =random.uniform(0, 1)
q = random.uniform(0, 1)
delta_x = r_min*np.power((r_max/r_min),p)*np.cos(2*np.pi*q)
delta_y = r_min*np.power((r_max/r_min),p)*np.cos(2*np.pi*q)
pos_q = np.array((pos[0]+delta_x, pos[1]+delta_y))

# Boundary conditions

if pos_q[0]>L:
	pos_q[0] = pos_q[0] - L
if pos_q[0]<0:
	pos_q[0] = pos_q[0] + L
if pos_q[1]>L:
	pos_q[1] = pos_q[1] - L
if pos_q[1]<0:
	pos_q[1] = pos_q[1] + L

