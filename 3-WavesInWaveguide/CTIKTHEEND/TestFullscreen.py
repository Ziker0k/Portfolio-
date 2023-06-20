#!/usr/bin/env python

# import usefull modules
import matplotlib
from numpy import *
from pylab import *
from scipy.integrate import ode
from scipy import interpolate
from scipy.integrate import ode

# Параметры волновода
m = 1
n = 0
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
time = 3 * (T / 40)
# Количество линий уровня
k = 10


def Formules(x, y, z, typeVolni, typePolya):
    if typeVolni == "TE":
        if typePolya == "E":
            Xcomp = -cos(m * pi * x / a) * sin(n * pi * y / b) * sin(w * time - beta * z)
            Ycomp = sin(m * pi * x / a) * cos(n * pi * y / b) * sin(w * time - beta * z)
            Zcomp = x * y * z * 0
            return Xcomp, Ycomp, Zcomp
        elif typePolya == "H":
            Xcomp = -sin(m * pi * x / a) * cos(n * pi * y / b) * sin(w * time - beta * z)
            Ycomp = -cos(m * pi * x / a) * sin(n * pi * y / b) * sin(w * time - beta * z)
            Zcomp = cos(m * pi * x / a) * cos(n * pi * y / b) * cos(w * time - beta * z)
            return Xcomp, Ycomp, Zcomp


xx = np.linspace(0, 23, 100)
yy = np.linspace(0, 10, 100)
X, Y = np.meshgrid(xx, yy)
Ex, Ey, Ez = Formules(X, Y, 0, "TE", "E")
Eabs = sqrt(Ex**2 + Ey**2)

# plt.figure(figsize=(4,4),facecolor="w")
# plt.streamplot(X, Y, Ex, Ey, color='r', density=1)
# plt.show()

# interpolate function of the Bx and Bz as functions of (x,z) position
fex = interpolate.interp2d(xx,yy,Ex)
fey = interpolate.interp2d(xx,yy,Ey)
def B_dir(t,p,fx,fy):
    ex = fx(p[0],p[1])
    ey = fy(p[0],p[1])
    n = sqrt(ex**2+ey**2)
    return [ex/n, ey/n]

print(B_dir)

# set the starting point of the magnetic field line
xstart = np.linspace(0, 23, 50)
ystart = np.linspace(0, 10, 50)
Xstart, Ystart = np.meshgrid(xstart, ystart)
places = np.vstack([Xstart,Ystart])

R=0.01
dt=0.8*R

# plot area
x0, x1=0, 23
y0, y1=0, 10

# set the ode function
r=ode(B_dir)
r.set_integrator('vode')
r.set_f_params(fex,fey)

xs,ys = [],[]
for p in places:
    x=[p[0]]
    y=[p[1]]
    r.set_initial_value([x[0], y[0]], 0)
    while r.successful():
        r.integrate(r.t+dt)
        x.append(r.y[0])
        y.append(r.y[1])
        hit_electrode=False
        # check if field line left drawing area
        if (not (x0<r.y[0] and r.y[0]<x1)) or (not (y0<r.y[1] and r.y[1]<y1)):
            break
    xs.append(x)
    ys.append(y)

plt.figure(figsize=(7,7),facecolor="w")

for x,y in zip(xs,ys):
    plt.plot(x,y,color="r")

plt.show()


