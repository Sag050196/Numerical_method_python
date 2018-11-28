"""

Here is the implementation of Finite Difference Time Domain in 1 dimension
Vanishing of the electric field in vacuum

"""
import numpy as np
import matplotlib.pyplot as plt
import math as m


def checker(e, n):
    ch = 0
    for i in range(n):  # Цикл бегающей по массиву значений эл поля
        if abs(e[i]) < 10**(-2):  # Условия, определяющие достигнутое значение эл поля
            ch = False
        else:
            return True
    return ch


N = 1000
sigma = 10.0
MaxTime = 500
delay = 70

#field
hy = np.zeros(N)
ez = np.zeros(N)

imp0 = 377.0#impedance

# mesh
x = np.linspace(0, N-1, N)

plt.ion()

for t in range(MaxTime):
    hy[:-1] += (ez[1:] - ez[:-1]) / imp0
    ez[1:] += (hy[1:] - hy[:-1]) * imp0
    # ez[int(N/2)-1] += m.exp(-((t-delay)/sigma)**2)

    if t == 1:  # Хевисайд в качестве источника
        ez[int(N/2)-1]=1
    if t < delay:
        ez[int(N / 2) - 1] += m.exp(-((t - delay) / sigma) ** 2)

    if t % 10:
        plt.plot(x, ez, color='red', linestyle='solid', linewidth='1', label='Ez')
        plt.plot(x, hy*imp0, color='blue', linestyle='dashed', label='Hy')
        plt.legend()
        #plt.savefig("img-%i.png" % t)
        plt.pause(0.00001)
        plt.clf()
        plt.ioff()
    if not checker(ez, N) and t > delay:
        plt.figure()
        plt.clf()
        plt.subplot(211)  # 2 строки, 1 столбец, первая строка
        plt.plot(x, ez, color='red', linestyle='solid', linewidth='1', label='Ez')
        plt.legend()
        plt.subplot(212)  # 2 строки, 1 столбец, вторая строка
        plt.plot(x, hy * imp0, color='blue', linestyle='dashed', label='Hy')
        plt.legend()
        plt.show()
        break
