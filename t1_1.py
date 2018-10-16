import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D


#Boundary function

def phix(x):
    return np.sin(x)/np.sin(a)

def phiy(x):
    return np.sinh(x)/np.sinh(b)

# Set maximum iteration
maxIter = 500

# Set size

a = 15
b = 7


m = 50
n = 50
# Set Dimension
lenX = m
lenY = n


# Boundary condition
Utop = phix
Ubottom = 0
Uleft = 0
Uright = phiy

# Initial guess of interior grid
Uguess = 0.5

# Set colour interpolation and colour map
colorinterpolation = 50
colourMap = plt.cm.jet
 # colourMap = plt.cm.coolwarm

# Set meshgrid
x = np.linspace(0, a, lenX)
y = np.linspace(0, b, lenY)

X, Y = np.meshgrid(x, y)
# Set array size and set the interior value with Tguess
U = np.empty((lenY, lenX))
U1 = np.empty((lenY, lenX))
U.fill(Uguess)

# Set Boundary condition
for i in range(lenX):
    U[(lenY-1):, i] = Utop(i*a/lenX)


U[:1, :] = Ubottom
for i in range(lenY):
    U[i, (lenX-1):] = Uright(i*b/lenY)

U[:, :1] = Uleft

# Iteration (We assume that the iteration is convergence in maxIter = 500)
print("Please wait for a moment")
for iteration in range(0, maxIter):
    for i in range(1, lenY-1):
        for j in range(1, lenX-1):
            U[i, j] = 1/2 *1/((lenX/a)**2 + (lenY/b)**2) * ((U[i+1][j] + U[i-1][j])/(b/lenY)**2 + (U[i][j+1] + U[i][j-1])/(a/lenX)**2)
print("Iteration finished")

print("Start solve analitic")
for i in range(0,lenY):
    for j in range(0,lenX):
        U1[i,j] = i*j/(a*b)*(a*b/(lenY*lenX))
        for n in range(1,55):
            U1[i,j] = U1[i,j]+ (2*((-1)**n) * (b**2)/(np.pi* n*((b**2) + (n**2) * (np.pi**2)) ) * np.sinh(math.pi*n*j*a/lenX/b)) * math.sin(np.pi * n* i*b/lenY/b)/np.sinh(np.pi*n*a/b) + (2*((-1)**n) * (a**2)/(np.pi* n*(a**2 - n**2 * np.pi**2) ) * np.sinh(np.pi*n*i*b/lenY/a)) * np.sin(np.pi * n* j*a/lenX/a)/np.sinh(np.pi*n*b/a)
print("Finished")

#fig = plt.figure()
#ax = Axes3D(fig)

# Configure the contour
plt.title("Solution explicit method")
plt.contourf(X, Y, U, colorinterpolation, cmap=colourMap)
#ax.plot_surface(X, Y, U, rstride=1, cstride=1)


# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()


#fig = plt.figure()
#ax = Axes3D(fig)
# Configure the contour
plt.title("Solution analitic ")
plt.contourf(X, Y, U1, colorinterpolation, cmap=colourMap)
#ax.plot_surface(X, Y, U1, rstride=1, cstride=1)
# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()

print("")
