import numpy as np
import matplotlib.pyplot as plt
import math
import importlib
import time
from mpl_toolkits.mplot3d import Axes3D



#Boundary function

def phix(x):
    return np.sin(x)/np.sin(a)

def phiy(x):
    return np.sinh(x)/np.sinh(b)

def SolvJacobi(a,m,b,n,maxIter):
    # Set Dimension
    lenX = m
    lenY = n


    # Boundary condition
    Utop = phix
    Ubottom = 0
    Uleft = 0
    Uright = phiy

    # Initial guess of interior grid
    Uguess = 0


    # Set meshgrid
    x = np.linspace(0, a, lenX)
    y = np.linspace(0, b, lenY)

    X, Y = np.meshgrid(x, y)
    # Set array size and set the interior value with Tguess
    U = np.empty((lenY, lenX))
    U.fill(Uguess)

    # Set Boundary condition

    U[(lenY-1), :] = Utop(x)
    U[:1, :] = Ubottom
    U[:, (lenX-1)] = Uright(y)
    U[:, :1] = Uleft
    for interation in range(0,maxIter):
       U[1:lenY-1,1:lenX-1] = 1/2 *(1/(((lenX-1)/a)**2 + ((lenY-1)/b)**2)) * ((U[2:,1:lenX-1] + U[:lenY-2,1:lenX-1])/(b/(lenY-1))**2 + (U[1:lenY-1,2:] + U[1:lenY-1,:lenX-2])/(a/(lenX-1))**2)
    return U

def SolvingNaA(a,m,b,n):
    # Set Dimension
    lenX = m
    lenY = n


    # Boundary condition
    Utop = phix
    Ubottom = 0
    Uleft = 0
    Uright = phiy

    #Convergence
    #tol = 2*10**(-4)

    # Initial guess of interior grid
    Uguess2 = 0
    # Set colour interpolation and colour map
    colorinterpolation = 50
    colourMap = plt.cm.jet
     # colourMap = plt.cm.coolwarm


    hx = a/(lenX-1)
    hy = b/(lenY-1)


        # Set meshgrid
    x = np.linspace(0, a, lenX)
    y = np.linspace(0, b, lenY)
    X, Y = np.meshgrid(x, y)
        # Set array size and set the interior value with Tguess
    Unew = np.empty((lenY, lenX))
    U1 = np.empty((lenY, lenX))
    Unew.fill(Uguess2)

        # Set Boundary condition



    Unew[lenY-1,:] = Utop(x)
    Unew[0, :] = Ubottom
    Unew[:, lenX-1] = Uright(y)
    Unew[:, 0] = Uleft
        #print("Please wait for a Analytic")
    for i in range(0,lenY):
        for j in range(0,lenX):
            U1[i,j] = (np.sin(j*a/(lenX-1)) * np.sinh(i*b/(lenY-1)))/(np.sin(a) * np.sinh(b))
                #print("Finished")
        # Iteration (We assume that the iteration is convergence in maxIter = 500)
        #print("Please wait for a moment")
        #while np.all(abs(U - Uold))>tol:



    alpha = (hy/hx)**2
    k = (n-2)*(m-2)
    A = np.zeros((k, k))
    N = m
    I = np.identity(N - 2)
    B = np.zeros((N - 2, N - 2))
    np.fill_diagonal(B, -2 * (1 + alpha))
    for i in range(N - 3):
        B[i + 1, i] = B[i, i + 1] = alpha
    for i in range(0, k, N - 2):
        A[i:i + N - 2, i:i + N - 2] = B

    for i in range(0, k - N+1, N - 2):
        A[i + N - 2:i + N - 4 + N, i:i + N - 2] = A[i:i + N - 2, i + N - 2:i + N - 4 + N] = I
    f = np.zeros((n, m))
    F = np.zeros(k)
    t = 0



    for i in range(1,n-1):
        for j in range(1,m-1):
            F[t] = (hx) ** 2 * f[i, j]
            if i-1 == 0:
                F[t] -= Unew[i-1, j]
            elif i+1 == n-1:
                  F[t] -= Unew[i + 1, j]
            if j+1 == m-1:
                F[t] -= alpha * Unew[i, j + 1]
            elif j-1 == 0:
                  F[t] -= alpha * Unew[i, j - 1]
            t += 1

    x = np.linalg.solve(A, F)
    Unew[1:n-1, 1:m-1] = x.reshape(n-2, m-2)
    return X,Y,Unew,U1



maxIter = 500
#Area of solutions
a = 15
b = 7
#Set of mersh
mx = 20
ny = 16

X,Y,Unew,U1 = SolvingNaA(a,mx,b,ny)
U = SolvJacobi(a,mx,b,ny,maxIter)



fig = plt.figure(1)
ax = Axes3D(fig)
plt.title("Numerical Solution  ")
#plt.contourf(X, Y, Unew, colorinterpolation, cmap=colourMap)
ax.plot_surface(X, Y, Unew, rstride=1, cstride=1)

# fig = plt.figure(2)
# ax = Axes3D(fig)
# plt.title("Numerical Solution Jacobi ")
# #plt.contourf(X, Y, Unew, colorinterpolation, cmap=colourMap)
# ax.plot_surface(X, Y, U, rstride=1, cstride=1)


fig = plt.figure(3)
ax = Axes3D(fig)
# Configure the contour
plt.title("Solution analitic ")
#plt.contourf(X, Y, U1, colorinterpolation, cmap=colourMap)
ax.plot_surface(X, Y, U1, rstride=1, cstride=1)
# Set Colorbar
#plt.colorbar()

fig = plt.figure(4)
ax = Axes3D(fig)
# Configure the contour
plt.title("ERROR for inverse matrix method")
#plt.contourf(X, Y, U, colorinterpolation, cmap=colourMap)
ax.plot_surface(X, Y, abs(U1-Unew), rstride=1, cstride=1)

# Show the result in the plot window
plt.show()


err1 = []
err2 = []
for mx in range(10,60):
    X,Y,Unew,U1 = SolvingNaA(a,mx,b,ny)
    U = SolvJacobi(a,mx,b,ny,maxIter)
    err1.append(np.sqrt(np.sum((U1 - Unew) ** 2) / mx / ny))
    err2.append(np.sqrt(np.sum((U1 - U) ** 2) / mx / ny))

plt.title("ERROR from iterations ")
plt.plot(range(10,60),err1,label = 'Revers matrix method')
plt.plot(range(10,60),err2,label = 'Iterational metod')
plt.legend()
plt.xlabel('Number of points on X')
plt.show()

print("")
