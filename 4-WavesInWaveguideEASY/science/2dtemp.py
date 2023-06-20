import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from numpy import (sin, cos, pi, sqrt, arange, meshgrid)
import seaborn as sns

# Параметры волновода
m = 1
n = 1
a = 23  # Размер волновода по x
b = 10  # Размер волновода по y
# Размер волновода по z

cc = 3 * (10 ** 11)  # Скорость света 3 мм/с
lambda_kr = (2 * a * b) / sqrt(pow(m * b, 2) + pow(n * a, 2))
lambda1 = lambda_kr - lambda_kr / 5
frequency = cc / lambda1
w = 2 * pi * frequency

beta = (w / cc) * (sqrt(1 - pow((lambda1 / lambda_kr), 2)))
kappa = sqrt((m * pi / a) ** 2 + (n * pi / b) ** 2)

# Отсчет
h = 1
# Время
T = 1 / frequency
time = 20 * T / 40
# Количество линий уровня
k = 10


def Formules(x, y, z, typeVolni, typePolya):
    if typeVolni == "TE":
        if typePolya == "E":
            Xcomp = -(w / (kappa ** 2)) * (n * pi / b) * cos(m * pi * x / a) * sin(
                n * pi * y / b) * sin(w * time - beta * z)
            Ycomp = (w / (kappa ** 2)) * (m * pi / a) * sin(m * pi * x / a) * cos(
                n * pi * y / b) * sin(w * time - beta * z)
            Zcomp = x * y * z * 0
            return Xcomp, Ycomp, Zcomp
        elif typePolya == "H":
            Xcomp = -(beta / (kappa ** 2)) * (m * pi / a) * sin(m * pi * x / a) * cos(
                n * pi * y / b) * sin(w * time - beta * z)
            Ycomp = -(beta / (kappa ** 2)) * (n * pi / b) * cos(m * pi * x / a) * sin(
                n * pi * y / b) * sin(w * time - beta * z)
            Zcomp = cos(m * pi * x / a) * cos(n * pi * y / b) * cos(w * time - beta * z)
            return Xcomp, Ycomp, Zcomp
    else:
        if typePolya == "E":
            Xcomp = (-(beta * m * pi) / ((kappa ** 2) * a)) * cos((m * pi * x) / a) * sin(
                (n * pi * y) / b) * cos(w * time - beta * z)
            Ycomp = (-(beta * n * pi) / ((kappa ** 2) * b)) * sin((m * pi * x) / a) * cos(
                (n * pi * y) / b) * cos(w * time - beta * z)
            Zcomp = sin((m * pi * x) / a) * sin((n * pi * y) / b) * cos(w * time - beta * z)
            return Xcomp, Ycomp, Zcomp
        elif typePolya == "H":
            Xcomp = -((w * n * pi) / ((kappa ** 2) * b)) * sin((m * pi * x) / a) * cos(
                (n * pi * y) / b) * sin(w * time - beta * z)
            Ycomp = ((w * m * pi) / ((kappa ** 2) * a)) * cos((m * pi * x) / a) * sin(
                (n * pi * y) / b) * sin(w * time - beta * z)
            Zcomp = 0 * x * y * z
            return Xcomp, Ycomp, Zcomp


def exp(z, time):
    return np.exp(w*time -1j * beta * z)

def IDiNahoyE(x, y, z):
    return np.real((w / (kappa ** 2 * 1j)) * (n * pi / b) * cos(m * pi * x / a) * sin(n * pi * y / b) * exp(z, time)),\
        np.real(-(w / (kappa ** 2 * 1j)) * (m * pi / a) * sin(m * pi * x / a) * cos(n * pi * y / b) * exp(z, time))

def IDiNahoyH(x, y, z):
    return np.real((beta/kappa**2 * 1j) * (m*pi/a) * sin((m * pi * x) / a) * cos((n * pi * y) / b) * exp(z, time)), \
        np.real((beta/kappa**2 * 1j) * (n*pi/b) * cos((m * pi * x) / a) * sin((n * pi * y) / b) * exp(z, time)), \
        np.real(cos((m * pi * x) / a) * cos((n * pi * y) / b) * exp(z, time))


cSt = lambda1
fig = plt.figure()
XY = fig.add_subplot(3, 1, 1)
YZ = fig.add_subplot(3, 1, 2)
XZ = fig.add_subplot(3, 1, 3)
x = np.linspace(0, a, 10)
y = np.linspace(0, b, 10)
z = np.linspace(0, cSt, 10)

X, Y, Z = np.meshgrid(x, y, z)

zplane = 0

X1, Y1 = np.meshgrid(x, y)
X2, Y2 = np.meshgrid(x, b/2)

EXmax, EYmax, EZmax = Formules(X, Y, Z, "TE", "E")


EU, EV, EW = Formules(X1, Y1, zplane, "TE", "E")
HU, HV, HW = Formules(X1, Y1, zplane, "TE", "H")

maxE = np.max((sqrt(EXmax ** 2 + EYmax ** 2 + EZmax ** 2)))
module = sqrt(EU ** 2 + EV ** 2 + EW ** 2) / (maxE)

cs = XY.contour(X1, Y1, module,  colors='r')
blya = []
blyay = []
for item in cs.collections:
   for i in item.get_paths():
        v = i.vertices
        x = v[:, 0]
        y = v[:, 1]
        blya.append(np.max(x))
        blyay.append(np.max(y))
        print(blyay)

XX, YY = np.meshgrid(blya, blyay)
EXX, EYY, EZZ = Formules(XX, YY, zplane, "TE", "E")
XY.quiver(XX, YY, EXX, EYY, pivot='tip', color='r', scale_units="xy")

plt.show()
