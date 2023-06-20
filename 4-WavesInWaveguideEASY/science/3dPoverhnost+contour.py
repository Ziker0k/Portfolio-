import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from numpy import (sin, cos, pi, sqrt, arange, meshgrid)
import seaborn as sns

# Параметры волновода
m = 1
n = 0
a = 23  # Размер волновода по x
b = 10  # Размер волновода по y
# Размер волновода по z

cc = 3 * (10 ** 11)  # Скорость света 3 мм/с
lambda_kr = (2 * a * b) / sqrt(pow(m * b, 2) + pow(n * a, 2))
lambda1 = a * 1.3
c = 100
frequency = cc / lambda1
print(frequency)
w = 2 * pi * frequency
print(lambda_kr, lambda1)

beta = (w / cc) * (sqrt(1 - pow((lambda1 / lambda_kr), 2)))

# Отсчет
h = 1
# Время
t = 0

T = lambda1/cc
print(T)

# Количество линий уровня
k = 10


def TE10_Hx(xconst, y, z):
    return (beta * a / pi) * sin(pi * xconst / a) * sin(w * t - beta * z)


def TE10_Hz(x, y, zconst):
    return cos(pi * x / a) * sin(w * t - beta * zconst)


def TE10_Hy(x, z):
    return ((beta * a / pi) * sin(pi * x / a) * sin(w * t - beta * z)) * (cos(pi * x / a) * sin(w * t - beta * z))


def TE10_HyKakEy(x, z):
    return (w * a / pi) * (sin((pi * x) / a)) * (cos((w * t) - (beta * z))) / (1e10 * 2 * a)


def TE10_HyPlusEyContour(x, z):
    return abs((beta * a / pi) * sin(pi * x / a) * sin(w * t - beta * z)) * \
        abs(cos(pi * x / a) * sin(w * t - beta * z)) \
        / ((w * a / pi) * (cos((pi * x) / a)) * (cos((w * (t + T/4)) - (beta * z))) / (1e10 * a / 10))


def TE10_Ey(x, z):
    return (w * a / pi) * (sin((pi * x) / a)) * (sin((w * t) - (beta * z))) / (1e10 * 2 * a)


x = arange(0, a+1, h)
y = arange(0, b+1, h)
z = arange(0, c, h)

x2, z2 = meshgrid(x, z)
y1, z1 = meshgrid(y, z)

Hy = TE10_Hy(x2, z2)
Ey = TE10_Ey(x2, z2)

print(np.amax(Ey))
print(np.amax(Hy))

ax = plt.axes(projection='3d')

X, Y, Z = meshgrid(x, y, z)


#ax.plot_surface(x2, z2, TE10_Hy(x2,z2), cmap='coolwarm', alpha=0.9)
#ax.plot_surface(y1, z1, TE10_Hx(a/4, y1, z1))
ax.scatter(TE10_Hx(X, Y, Z))
#ax.contour(x2, Ey, z2, c, cmap='plasma', alpha=0.5)
#ax.set_xlim(0, c / 2)
#ax.set_ylim(0, c)
#ax.contour(x2, z2, TE10_HyKakEy(x2, z2), 5)
# ax.plot_surface(x2, z2, TE10_Ey(x2, z2), alpha=0.1, cmap='OrRd')
