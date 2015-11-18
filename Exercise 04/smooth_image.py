import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


h=10
def w(x):
    if x/h<1/2 and x/h>=0:
        l=1.0-6.0*np.power(x/h,2)+6.0*np.power(x/h,3)
    elif x/h>=1/2 and x/h<1:
        l=2.0*np.power((1-(x/h)),3)
    else:
        l=0
    return l 

# Normalize function

def integral(x):
    return w(x)*2*np.pi*x
result, err =quad(integral,0,h)

k=1/result



#Reads a square image in 8-bit/color PPM format from the given file. Note: No checks on valid format are done.
def readImage(filename):
    f = file(filename,"rb")
    
    f.readline()
    s = f.readline()
    f.readline()
    (pixel, pixel) = [t(s) for t,s in zip((int,int),s.split())]
    
    data = np.fromfile(f,dtype=np.uint8,count = pixel*pixel*3)
    img = data.reshape((pixel,pixel,3)).astype(np.double)
    
    f.close()
    
    return img, pixel
    

#Writes a square image in 8-bit/color PPM format.
def writeImage(filename, image):
    f = file(filename,"wb")
    
    pixel = image.shape[0]
    f.writelines("P6\n%d %d\n%d\n"%(pixel, pixel, 255))
    
    image = image.astype(np.uint8)
    
    image.tofile(f)
    
    f.close()
    
    
    
    

img, pixel = readImage("aq-original.ppm")


#Now we set up our desired smoothing kernel. We'll use complex number for it even though it is real. 
kernel_real = np.zeros((pixel,pixel),dtype=np.complex)

hsml = 10.

#start = [pixel-70,pixel-40]
start = [0,0]
counter =0
#now set the values of the kernel 
for i in np.arange(pixel):
    for j in np.arange(pixel):
        
        #TODO: do something sensible here to set the real part of the kernel
        r= np.sqrt(np.power(i-start[0],2)+np.power(j-start[0],2))
        kernel_real[i, j] = k*w(r) 
        if kernel_real[i, j]!=0:
            print kernel_real[i, j]
            counter = counter +1
#print counter

#Let's calculate the Fourier transform of the kernel
kernel_kspace = np.fft.fft2(kernel_real)


#further space allocations for image transforms
color_real = np.zeros((pixel,pixel),dtype=np.complex)


#we now convolve each color channel with the kernel using FFTs
for colindex in np.arange(3):
    #copy input color into complex array
    color_real[:,:].real = img[:,:,colindex]
    
    
    #forward transform
    color_kspace = np.fft.fft2(color_real)
    
    #multiply with kernel in Fourier space
    #TODO: fill in code here
    #color_kspace = np.real(color_kspace * kernel_kspace) #take only real part
    color_kspace = np.multiply(color_kspace,kernel_kspace)
    #backward transform
    color_real = np.fft.ifft2(color_kspace)
    
    
    #copy real value of complex result back into color array
    img[:,:,colindex] = color_real.real
    
writeImage("aq-smoothed.ppm", img)
