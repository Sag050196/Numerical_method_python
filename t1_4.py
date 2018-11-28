from numpy import *
import matplotlib.pyplot as plt

def ExactSolution(x,t):
  return outer(sin(pi*x),1/pi*sin(pi*t) + cos(pi*t))

def Solver(x,t):
  #initialization
  tau = t[1] - t[0]
  h = x[1] - x[0]
  Nx = len(x)
  Nt = len(t)
  U = zeros((Nx,Nt))

  #initial condition
  U[:,0] = sin(pi*x)
  U[:,1] = tau*sin(pi*x) + U[:,0]

  #matrix
  c = tau**2/h**2
  print(c)
  A = 2*(1-c)*identity(Nx-2) + c*eye(Nx-2,Nx-2,1) + c*eye(Nx-2,Nx-2,-1)
  for k in range(1,Nt-1):
    U[1:-1,k+1] = A.dot(U[1:-1,k]) - U[1:-1,k-1]
  return U




x = linspace(0, 1, 32)
t = linspace(0, 10, 310)
Uexplicit = Solver(x, t)
exactU = ExactSolution(x,t)
plt.pcolormesh(t, x, exactU )
plt.title('Analitic solution')
plt.colorbar()
plt.show()

plt.pcolormesh(t, x, Uexplicit )
plt.title('Method solution')
plt.colorbar()
plt.show()
# for i in range(0, k+1,int(k/k)):
#    plt.plot(x, u[i, :])
#    plt.show()
plt.pcolormesh(t, x, abs(exactU-Uexplicit) )

plt.title('Error')
plt.colorbar()
plt.show()


#Main
tmax = 1
Nx = 10
N = arange(5,1000,10)
TT = []
errExplicit = []
Courant = []
for Nt in N:
  x = linspace(0,1,Nx)
  t = linspace(0,tmax,Nt)
  tau = t[1] - t[0]
  h = x[1] - x[0]
  TT.append(tau)
  #print(tau**2 / h**2)
  Courant.append(tau**2/h**2)
  Uexplicit = Solver(x, t)
  exactU = ExactSolution(x,t)
  errExplicit.append(sqrt(sum((Uexplicit - exactU)**2)/Nx/Nt))

#ploting
plt.loglog(Courant,errExplicit)
plt.ylabel('error')
plt.xlabel('Courant number')
plt.show()
