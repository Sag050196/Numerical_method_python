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
        if abs(e[i]) < 10**(-4):  # Условия, определяющие достигнутое значение эл поля
            ch = False
        else:
            return True
    return ch


N = 100
sigma = 10.0
MaxTime = 500
delay = 50


n=2;
Courant = 1*n;

hy = np.zeros(N)
ez = np.zeros(N)

imp0 = 377.0;

# Иницилизация координат пространства для построения полей
x = np.linspace(0, N-1, N)

plt.ion()

for t in range(MaxTime):
    hy[-1] = hy[-2]
    hy[:-1] += (ez[1:] - ez[:-1]) / imp0 * Courant
    ez[0] = ez[1]
    ez[1:] += (hy[1:] - hy[:-1]) * imp0 * Courant /n**2
    ez[int(N/2)-1] += 1/n**2 *  m.exp(-((t-delay)/sigma)**2)

    # if t == 1:  # Хевисайд в качестве источника
    #      ez[int(N/2)-1]=1
    # if t < delay:
    #     ez[int(N / 2) - 1] += m.exp(-((t - delay) / sigma) ** 2)

    if t % 20:
        plt.plot(x, ez, color='red', linestyle='solid', linewidth='1', label='Ez')
        plt.plot(x, hy*imp0, color='blue', linestyle='dashed', label='Hy')
        plt.legend()
        plt.pause(0.001)
        plt.clf()

    plt.ioff()

    if not checker(ez, N) and t > delay:
        plt.clf()
        plt.subplot(211)  # 2 строки, 1 столбец, первая строка
        plt.plot(x, ez, color='red', linestyle='solid', linewidth='1', label='Ez')
        plt.subplot(212)  # 2 строки, 1 столбец, вторая строка
        plt.plot(x, hy * imp0, color='blue', linestyle='dashed', label='Hy')
        plt.show()
        break  # Прерывание цикла и, как следдствие, завершение программы
