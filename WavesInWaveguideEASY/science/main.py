import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from numpy import (sin, cos, pi, sqrt, arange, meshgrid)

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
t = 6
# Количество линий уровня
k = 10


# Алгоритм определения констант для проекций (где будет сделан срез)


# # H10
# def TE10_H_XY(x, y):
#     return abs(sin(kappaX * x)) * exp(hh * tan(w * t) * y)
#
#
# def TE10_H_XZ(x, z):
#     return abs(sin(kappaX * x)) * cos(w * t - hh * z)
#
#
# def TE10_H_YZ(y, z):
#     return exp(kappaX * (cos(kappaX * C1) / sin(kappaX * C1)) * y) * cos(w * t - hh * z)
#
#
# # E10
# def TE10_E_XY(x, y):
#     return (w / (kappaX * cc)) * abs(sin(kappaX * x)) * sin(w * t)

def makeData1(x, y, z):
    X = arange(0, x, h)
    Y = arange(2*y, 3*y, h)
    Z = arange(0, z, h)
    xgrid, ygrid, zgrid = meshgrid(X, Y, Z)

    return xgrid, ygrid, zgrid

def makeData(x, y, z):
    X = arange(0, x, h)
    Y = arange(0, y, h)
    Z = arange(0, z, h)
    xgrid, ygrid, zgrid = meshgrid(X, Y, Z)
    return xgrid, ygrid, zgrid


def makeDataWithTwo(a1, a2):
    b1 = arange(0, a1, h)
    b2 = arange(0, a2, h)
    a1grid, a2grid = meshgrid(b1, b2)
    return a1grid, a2grid


fig = plt.figure()


def plotting():
    def TE10_Ey(x, y, z):
        # return (w / (kappaX * cc)) * abs(sin(kappaX * x)) * sin(w * t)
        return abs((w * a / pi)) * abs(sin((pi * x) / a)) * abs(sin((w * t) - (beta * z))) * b / (1e10 * 2 * a)

    def TE10_Hx(x, y, z):
        return ((beta*a/pi) * sin(pi*x/a) * sin(w*t - beta*z))

    def TE10_Hz(x, y, z):
        return cos(pi*x/a) * sin(w*t - beta*z)

    def TE10_Hy(x, y, z):
        return (cos(pi*x/a) * sin(w*t - beta*z)) * ((beta*a/pi) * sin(pi*x/a) * sin(w*t - beta*z))

    def gran(x, z):
        return b

    #x2, z2 = makeData(a, c)
    #
    # Ey = TE10_Ey(x2, z2)
    # ax = plt.axes(projection='3d')
    #
    # ax.plot_surface(x2, z2, Ey, cmap='plasma')
    # ax.scatter(x2, z2, gran(x2, z2), alpha=0)
    #
    # xLabel = ax.set_xlabel('\nXXX xxxxxx xxxx x xx x')
    # yLabel = ax.set_ylabel('\nZZZzzzzz')
    # zLabel = ax.set_zlabel('\nYyyyy yyyy y')
    # ax.set_xlim(0, c)
    # ax.set_ylim(0, c)
    # plt.show()


    # ax1 = plt.axes(projection='3d')
    # y,z = makeData(b, c)
    # x1,y1 = makeData(a, b)
    #
    # print(Hy)
    # ax1.plot_surface(x2, z2, Hy, cmap='inferno')

    #
    # ax1.set_xlabel('\nXXX xxxxxx xxxx x xx x')
    # ax1.set_ylabel('\nZ')
    # ax1.set_zlabel('\nY')

    # import seaborn

    # load "flights" dataset
    #data = np.random.random((12, 12))
    # creating a dummy dataset


    ax = plt.axes(projection='3d')
    x, y, z = makeData(a, b, c)
    x1, y1, z1 = makeData1(a, b, c)

    Value = TE10_Ey(x, y, z)
    Value1 = TE10_Hx(y1, x1, z1)
    Value2 = TE10_Hy(x1, y1, z1)

    granicax = [0, c/2]
    granicay = [0, c/2]
    granicaz = [0, c]


    #ax.scatter(x, y, z, c=Value, lw=0, s=20)
    ax.scatter(x1, y1, z1, c=Value2, lw=0, s=10, cmap='Blues', marker='.')
    ax.scatter(granicax, granicay, granicaz, alpha=0)

    ax.set_xlabel('\nX')
    ax.set_ylabel('\nY')
    ax.set_zlabel('\nZ')

    color_map = matplotlib.cm.ScalarMappable()
    color_map.set_array(Value)
    plt.colorbar(color_map)
    plt.show()

plotting()
