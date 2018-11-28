import numpy as np
import matplotlib.pyplot as plt


def Solving(t_max,k,n):
    h = 1/(n)  #length of one space step
    tau = t_max/(k)  #length of one time step
    r = tau/h**2 #Stable factor should be less 1/2
    # print("Stable factor:", r)
    x = np.linspace(0, 1, n+1)
    ti = np.linspace(0, t_max, k+1)
    #initial condition
    U = np.zeros((k+1, n+1))
    U[0, :] = np.sin(np.pi*x)
    U[:, 0] = 0
    U[:, n] = 0
    Uext = np.zeros((k+1, n+1))
    for i in range(1, k):
        for j in range(1, n):
            Uext[i,j] = np.exp(-np.pi**2*ti[i])*np.sin(np.pi*x[j])
        #plt.plot(x,u_th)
        #plt.show()
    # scheme of solving
    # theta = 0 (explicit scheme), theta = 1 (implicit scheme), and theta = 1/2 (Crank-Nicolson scheme)
    theta = 1/2
    p = (n-1)
    #A = np.zeros((p, p))
    N = n
    B = np.zeros((N - 1, N - 1))
    np.fill_diagonal(B, -(1/tau+2*theta/h**2))
    for i in range(N - 2):
        B[i, i + 1]=B[i+1, i ] = theta/h**2
       # for i in range(0, p, N - 1):
            #A[i:i + N - 1, i:i + N - 1] = B

    #f = np.zeros((k+1, n+1))
    F = np.zeros(p)
    t = 0
    for i in range(1, k):
        F = np.zeros(p)
        t=0
        for j in range(1, n):
            F[t] = -(1/tau+2*(theta-1)/h**2)*U[i - 1, j]+(theta-1)/h**2*(U[i - 1,j-1]+U[i - 1, j+1])
            t += 1
        solve = np.linalg.solve(B, F)
        U[i, 1:n] = solve.reshape(1, n-1)
    return x,ti,U,Uext


k = 1000
n = 20
t_max = 1 #sec

x,ti,U,Uext = Solving(t_max,k,n)



plt.pcolormesh(x, ti, Uext )
plt.title('Analitic solution')
plt.colorbar()
plt.show()

plt.pcolormesh(x, ti, U )
plt.title('Method solution')
plt.colorbar()
plt.show()
# for i in range(0, k+1,int(k/k)):
#    plt.plot(x, u[i, :])
#    plt.show()
plt.pcolormesh(x, ti, abs(U-Uext) )

plt.title('Error')
plt.colorbar()
plt.show()

err = []
for k in range(100,1000,10):
    x,ti,U,Uext = Solving(t_max,k,n)
    err.append(np.sqrt(np.sum((U - Uext) ** 2) / k / n))

plt.title("ERROR from iterations ")
plt.loglog(range(100,1000,10),err)
plt.xlabel('Number of points')
plt.show()
