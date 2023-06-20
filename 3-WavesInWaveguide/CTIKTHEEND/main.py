import matplotlib
import tkinter.messagebox
import customtkinter
from customtkinter import *
from customtkinter import CTk
from tkinter import messagebox
import numpy as np
from numpy import (sin, cos, tan, pi, exp, sqrt, arange, meshgrid, linspace)
from threading import Thread
import time as tm
import re

matplotlib.use('TkAgg')
from matplotlib import pyplot as plt, cm, pyplot
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.pyplot import (clf, gcf, subplots_adjust, figure)

customtkinter.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("dark-blue")  # Themes: "blue" (standard), "green", "dark-blue"


def testVal(inStr, acttyp):
    if acttyp == '1':  # insert
        if not inStr.isdigit():
            return False
    return True


def get(entry):
    value = entry.get()
    try:
        return float(value)
    except ValueError:
        return None


def click():
    if (len(root.moda_m_entry.get()) == 0) or (len(root.moda_n_entry.get()) == 0) or (
            len(root.size_a_entry.get()) == 0) or (len(root.size_b_entry.get()) == 0) \
            or (len(root.size_b_entry.get()) == 0) or (len(root.frequency_entry.get()) == 0):
        tkinter.messagebox.showerror(title=None, message="Заполните все поля!")
        return

    if root.type_volni_entry.get() == "TM" and (root.moda_m_entry.get() == '0' or root.moda_n_entry.get() == '0'):
        tkinter.messagebox.showerror(title=None, message="Мода ТМ может существовать только при m>0 и n>0!")
        return

    cc = 3 * (10 ** 11)  # Скорость света 3 мм/с

    # Параметры волновода
    m = get(root.moda_m_entry)
    n = get(root.moda_n_entry)
    a = get(root.size_a_entry)  # Размер волновода по x в мм
    b = get(root.size_b_entry)  # Размер волновода по y в мм
    c = get(root.size_c_entry)  # Размер волновода по z в мм

    frequency_ggz = get(root.frequency_entry)  # Частота в ГГц
    frequency = frequency_ggz * 1e9

    lambda_kr = (2 * a * b) / sqrt(pow(m * b, 2) + pow(n * a, 2))
    lambda1 = cc / frequency
    w = 2 * pi * frequency

    f_kr = cc / lambda_kr

    if frequency < f_kr:
        tkinter.messagebox.showerror(title=None, message="Частота ниже критической!")
        root.label.configure(text="Критическая частота(ГГц): " + str(f_kr / 1e9))
        return

    root.label.configure(text="Критическая частота(ГГц): " + str(f_kr / 1e9))

    if root.type_polya_entry.get() == "Модуль":
        window = Window()
    else:
        window = Vector()


root = CTk()
root.title("Волновод")
root.geometry(f"{600}x{600}")
root.update_idletasks()
w, h, sx, sy = map(int, re.split('x|\+', root.winfo_geometry()))
sw = (root.winfo_rootx() - sx) * 2 + w
sh = (root.winfo_rooty() - sy) + (root.winfo_rootx() - sx) + h
sx = (root.winfo_screenwidth() - sw) / 2
sy = (root.winfo_screenheight() - sh) / 2
root.wm_geometry('+%d+%d' % (sx, sy))

# configure grid layout (4x4)

root.grid_columnconfigure(0, weight=1)
root.grid_rowconfigure((0, 1, 2, 3), weight=1)

# create sidebar frame with widgets
root.sidebar_frame = customtkinter.CTkFrame(root, width=300, corner_radius=0)
root.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
root.sidebar_frame.grid_rowconfigure((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), weight=1)
root.sidebar_frame.grid_columnconfigure((0, 1), weight=1)

root.logo_label = customtkinter.CTkLabel(root.sidebar_frame, text="Характеристики волновода",
                                         font=customtkinter.CTkFont(size=20, weight="bold"))
root.logo_label.grid(row=0, columnspan=2, padx=20, pady=(20, 10))

# row1
root.label_n = customtkinter.CTkLabel(root.sidebar_frame, text="Выберите тип волны:")
root.label_n.grid(row=1, column=0, sticky="w", padx=15)
root.type_volni_entry = customtkinter.CTkOptionMenu(root.sidebar_frame, values=["TE", "TM"], fg_color="#565B5E",
                                                    button_color="#565B5E", button_hover_color="#464A4D",
                                                    dropdown_fg_color="#464A4D", dropdown_hover_color="#464A4D")
root.type_volni_entry.grid(row=1, column=1, padx=15)

# row2
root.moda_label = customtkinter.CTkLabel(root.sidebar_frame, text="Выберите моду:")
root.moda_label.grid(row=2, column=0, sticky="w", padx=15, pady=(15, 0))

# row3
root.moda_m = customtkinter.CTkLabel(root.sidebar_frame, text="m: ")
root.moda_m.grid(row=3, column=0, sticky="w", padx=15)
root.moda_m_entry = customtkinter.CTkEntry(root.sidebar_frame, width=50)
root.moda_m_entry['validatecommand'] = (root.moda_m_entry.register(testVal), '%P', '%d')
root.moda_m_entry.grid(row=3, column=0, padx=(30, 0), sticky="w")
root.moda_n = customtkinter.CTkLabel(root.sidebar_frame, text="n: ")
root.moda_n.grid(row=3, column=1, sticky="w", padx=15)
root.moda_n_entry = customtkinter.CTkEntry(root.sidebar_frame, width=50)
root.moda_n_entry['validatecommand'] = (root.moda_n_entry.register(testVal), '%P', '%d')
root.moda_n_entry.grid(row=3, column=1, padx=(30, 0), sticky="w")

# row4
root.frequency_label = customtkinter.CTkLabel(root.sidebar_frame, text="Введите частоту(ГГц): ")
root.frequency_label.grid(row=4, column=0, padx=15, pady=15, sticky="w")
root.frequency_entry = customtkinter.CTkEntry(root.sidebar_frame)
root.frequency_entry['validatecommand'] = (root.frequency_entry.register(testVal), '%P', '%d')
root.frequency_entry.grid(row=4, column=1, sticky="w", padx=15)

root.label = customtkinter.CTkLabel(root.sidebar_frame, text="")
root.label.grid(row=5, columnspan=2, padx=15, pady=15, sticky="SW")

# row6
root.size_label = customtkinter.CTkLabel(root.sidebar_frame, text="Введите размеры волновода(мм):")
root.size_label.grid(row=6, column=0, sticky="w", padx=15)

# row7
root.size_frame = customtkinter.CTkFrame(root.sidebar_frame, corner_radius=0)
root.size_frame.grid(row=7, columnspan=2, sticky="nsew")
root.size_frame.grid_columnconfigure((0, 1, 2), weight=1)
root.size_frame.grid_rowconfigure((0, 1), weight=1)

root.size_a_label = customtkinter.CTkLabel(root.size_frame, text="a:")
root.size_a_label.grid(row=0, column=0, sticky="nsew", padx=15)
root.size_a_entry = customtkinter.CTkEntry(root.size_frame, width=50)
root.size_a_entry['validatecommand'] = (root.size_a_entry.register(testVal), '%P', '%d')
root.size_a_entry.grid(row=1, column=0, padx=15, pady=(0, 15), sticky="nsew")

root.size_b_label = customtkinter.CTkLabel(root.size_frame, text="b:")
root.size_b_label.grid(row=0, column=1, sticky="nsew", padx=15)
root.size_b_entry = customtkinter.CTkEntry(root.size_frame, width=50)
root.size_b_entry['validatecommand'] = (root.size_a_entry.register(testVal), '%P', '%d')
root.size_b_entry.grid(row=1, column=1, padx=15, pady=(0, 15), sticky="nsew")

root.size_c_label = customtkinter.CTkLabel(root.size_frame, text="c:")
root.size_c_label.grid(row=0, column=2, sticky="nsew", padx=15)
root.size_c_entry = customtkinter.CTkEntry(root.size_frame, width=50)
root.size_c_entry['validatecommand'] = (root.size_a_entry.register(testVal), '%P', '%d')
root.size_c_entry.grid(row=1, column=2, padx=15, pady=(0, 15), sticky="nsew")

#####row8
root.type_polya_label = customtkinter.CTkLabel(root.sidebar_frame, text="Выберите тип отображения поля:")
root.type_polya_label.grid(row=8, column=0, padx=15)
root.type_polya_entry = customtkinter.CTkOptionMenu(root.sidebar_frame, values=["Модуль", "Вектор"], fg_color="#565B5E",
                                                    button_color="#565B5E", button_hover_color="#464A4D",
                                                    dropdown_fg_color="#464A4D", dropdown_hover_color="#464A4D")
root.type_polya_entry.grid(row=8, column=1, padx=15)

root.postroit = customtkinter.CTkButton(root.sidebar_frame, text="Построить", command=click,
                                        fg_color="#565B5E", hover_color="#464A4D") \
    .grid(row=9, columnspan=2)


class Window(CTk):
    time = 0
    check = False

    def __init__(self):

        super().__init__()

        # Параметры волновода
        m = get(root.moda_m_entry)
        n = get(root.moda_n_entry)
        a = get(root.size_a_entry)  # Размер волновода по x
        b = get(root.size_b_entry)  # Размер волновода по y
        frequencyGGz = get(root.frequency_entry)
        type_wave = root.type_volni_entry.get()
        self.type_polya = ""

        self.Xmesh = []
        self.Ymesh = []
        self.Zmesh = []
        self.kolvo_tochek = 0

        cc = 3 * (10 ** 11)  # Скорость света 3 мм/с
        lambda_kr = (2 * a * b) / sqrt(pow(m * b, 2) + pow(n * a, 2))
        lambda1 = cc / (frequencyGGz * 1e9)

        # c = lambda1 * (m + n)  # Размер волновода по z
        c = get(root.size_c_entry)  # Размер волновода по z

        frequency = (frequencyGGz * 1e9)
        w = 2 * pi * frequency
        beta = (w / cc) * (sqrt(1 - pow((lambda1 / lambda_kr), 2)))
        kappa = sqrt((m * pi / a) ** 2 + (n * pi / b) ** 2)

        f_kr = cc / lambda1
        T = 1 / f_kr

        Hm = 20
        Emax = 10000
        Eps = 1
        mu = 1

        def play():
            self.check = True
            while self.check:
                self.time += T / 40
                drawImg()
                tm.sleep(0.01)

        def time_minus():
            self.time -= T / 40
            drawImg()

        def time_plus():
            self.time += T / 40
            drawImg()

        def startPlay():
            if not self.check:
                thread = Thread(target=play, daemon=True)
                thread.start()

        def stopPlay():
            self.check = False

        def makeData(x, y, z, shag):
            X = np.linspace(0, x, int(x / shag))
            Y = np.linspace(0, y, int(y / shag))
            Z = np.linspace(0, z, int(z / shag))
            xgrid, ygrid, zgrid = meshgrid(X, Y, Z)

            return xgrid, ygrid, zgrid

        def calculation():
            cax.clear()
            k = get(k_entry)
            if get(k_entry) == 0:
                tkinter.messagebox.showerror(title=None, message="Заполните поле точек!")
                return
            self.type_polya = chose_field_entry.get()
            self.kolvo_tochek = int(get(k_entry))

            self.Xmesh, self.Ymesh, self.Zmesh = makeData(a, b, c, min(a, b, c) / self.kolvo_tochek)

            U, V, W = Formules(self.Xmesh, self.Ymesh, self.Zmesh, type_wave, self.type_polya)
            UVWmodule = (abs(U + V + W))

            mappable = cm.ScalarMappable(cmap='turbo')
            mappable.set_array(UVWmodule)
            self.cbar = fig.colorbar(mappable, cax=cax, orientation='vertical')

            drawImg()

        def Formules(x, y, z, typeVolni, typePolya):
            if typeVolni == "TE":
                if typePolya == "E":
                    Xcomp = -(w * mu / (kappa ** 2)) * (n * pi / b) * Hm * cos(m * pi * x / a) * sin(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Ycomp = ((w * mu) / (kappa ** 2)) * (m * pi / a) * Hm * sin(m * pi * x / a) * cos(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Zcomp = x * y * z * 0
                    return Xcomp, Ycomp, Zcomp
                elif typePolya == "H":
                    Xcomp = -(beta / (kappa ** 2)) * (m * pi / a) * Hm * sin(m * pi * x / a) * cos(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Ycomp = -(beta / (kappa ** 2)) * (n * pi / b) * Hm * cos(m * pi * x / a) * sin(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Zcomp = Hm * cos(m * pi * x / a) * cos(n * pi * y / b) * cos(w * self.time - beta * z)
                    return Xcomp, Ycomp, Zcomp
            else:
                if typePolya == "E":
                    Xcomp = (-(beta * m * pi) / ((kappa ** 2) * a)) * Emax * cos((m * pi * x) / a) * sin(
                        (n * pi * y) / b) * cos(w * self.time - beta * z)
                    Ycomp = (-(beta * n * pi) / ((kappa ** 2) * b)) * Emax * sin((m * pi * x) / a) * cos(
                        (n * pi * y) / b) * cos(w * self.time - beta * z)
                    Zcomp = Emax * sin((m * pi * x) / a) * sin((n * pi * y) / b) * cos(w * self.time - beta * z)
                    return Xcomp, Ycomp, Zcomp
                elif typePolya == "H":
                    Xcomp = -((w * Eps * n * pi) / ((kappa ** 2) * b)) * Emax * sin((m * pi * x) / a) * cos(
                        (n * pi * y) / b) * sin(w * self.time - beta * z)
                    Ycomp = ((w * Eps * m * pi) / ((kappa ** 2) * a)) * Emax * cos((m * pi * x) / a) * sin(
                        (n * pi * y) / b) * sin(w * self.time - beta * z)
                    Zcomp = 0 * x * y * z
                    return Xcomp, Ycomp, Zcomp

        def GetVolnovod(a, b, c):
            pogr = min(a, b, c) / 15
            xx = [0 - pogr, a + pogr]
            yy = [0 - pogr, b + pogr]
            zz = [0 - pogr, c + pogr]

            X = np.array(
                [xx[0], xx[1], xx[1], xx[0], xx[0], xx[0], xx[1], xx[1], xx[1], xx[1], xx[1], xx[1], xx[0], xx[0],
                 xx[0], xx[0]])
            Y = np.array(
                [yy[0], yy[0], yy[0], yy[0], yy[0], yy[1], yy[1], yy[0], yy[1], yy[1], yy[0], yy[1], yy[1], yy[0],
                 yy[1], yy[1]])
            Z = np.array(
                [zz[0], zz[0], zz[1], zz[1], zz[0], zz[0], zz[0], zz[0], zz[0], zz[1], zz[1], zz[1], zz[1], zz[1],
                 zz[1], zz[0]])
            return X, Y, Z

        # Dynamic drawing function
        def drawImg():
            if len(self.type_polya) != 0:
                self.tochki.remove()

                X, Y, Z = makeData(a, b, c, min(a, b, c) / self.kolvo_tochek)

                U, V, W = Formules(X, Y, Z, type_wave, self.type_polya)

                self.tochki = ax.scatter(X, Z, Y, c=(sqrt(U ** 2 + V ** 2 + W ** 2)),
                                         lw=0, s=1e3 * (m + n) / self.kolvo_tochek ** 2,
                                         marker='s', cmap='turbo', alpha=1)

                canvas.draw()

        # конфигурация окна
        self.title("График")
        self.w, self.h = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry("%dx%d" % (self.w, self.h))

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        frameOfGraph = customtkinter.CTkFrame(self, fg_color="lightgray", corner_radius=0)
        frameOfGraph.grid(row=0, column=0, sticky='nsew')

        frameOfGraph.grid_rowconfigure((0, 1, 2, 3), weight=2)
        frameOfGraph.grid_rowconfigure(4, weight=1)
        frameOfGraph.grid_columnconfigure((0, 1), weight=1)
        frameOfGraph.grid_columnconfigure(2, weight=3)

        chose_field_label = customtkinter \
            .CTkLabel(frameOfGraph, text="Выберите поле:", bg_color="transparent",
                      text_color="blue") \
            .grid(row=0, column=0)
        chose_field_entry = customtkinter.CTkOptionMenu(frameOfGraph, values=["E", "H"], fg_color="#565B5E",
                                                        button_color="#565B5E", button_hover_color="#464A4D",
                                                        dropdown_fg_color="#464A4D", dropdown_hover_color="#464A4D")
        chose_field_entry.grid(row=0, column=1)

        k_label = customtkinter \
            .CTkLabel(frameOfGraph, text="Введите количество точек\n(Влияет на скорость построения!) ",
                      bg_color="transparent", text_color="red") \
            .grid(row=1, column=0)
        k_entry = customtkinter.CTkEntry(frameOfGraph, width=50)
        k_entry['validatecommand'] = (k_entry.register(testVal), '%P', '%d')
        k_entry.grid(row=1, column=1)

        frameOfChanging = customtkinter.CTkFrame(frameOfGraph, fg_color="gray", corner_radius=0)
        frameOfChanging.grid(row=2, column=0, rowspan=3, columnspan=2, sticky='nsew')
        frameOfChanging.grid_rowconfigure((0, 1, 2), weight=1)
        frameOfChanging.grid_columnconfigure((0, 1), weight=1)

        CTkButton(frameOfChanging, text="Рассчитать", fg_color="#565B5E", hover_color="#464A4D", command=calculation) \
            .grid(row=0, columnspan=2)
        CTkButton(frameOfChanging, text="-", fg_color="#565B5E", hover_color="#464A4D", command=time_minus) \
            .grid(row=1, column=0)
        CTkButton(frameOfChanging, text='+', fg_color="#565B5E", hover_color="#464A4D", command=time_plus) \
            .grid(row=1, column=1)
        CTkButton(frameOfChanging, text="Play", fg_color="#565B5E", hover_color="#464A4D", command=startPlay) \
            .grid(row=2, column=0)
        CTkButton(frameOfChanging, text="Stop", fg_color="#565B5E", hover_color="#464A4D", command=stopPlay) \
            .grid(row=2, column=1)

        fig = figure(facecolor="white")
        ax = fig.add_subplot(projection='3d')
        g, h, j = GetVolnovod(a, b, c)
        ax.plot(g, j, h, color='r')
        cax = fig.add_subplot([0.03, 0.5, 0.02, 0.47])

        ax.set_xlim(0, max(a, b, c))
        ax.set_ylim(0, max(a, b, c))
        ax.set_zlim(0, max(a, b, c))

        self.tochki = ax.scatter(0, 0, 0, alpha=0)

        mappable = cm.ScalarMappable(cmap='turbo')
        self.cbar = fig.colorbar(mappable, cax=cax, orientation='vertical')

        canvas = FigureCanvasTkAgg(fig, master=frameOfGraph)
        canvas.get_tk_widget().grid(row=0, rowspan=4, column=2, sticky="nsew")

        toolbarFrame = CTkFrame(master=frameOfGraph, corner_radius=0, fg_color='white')
        toolbarFrame.grid(row=4, column=2, sticky='nsew')
        toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

        self.protocol("WM_DELETE_WINDOW", self.quit_clicked)
        self.mainloop()

    def quit_clicked(self):
        if messagebox.askokcancel("Закрыть", "Действительно хотите закрыть окно?"):
            self.check = False
            self.destroy()


class Vector(CTk):
    time = 0
    check = False
    visible = False

    def __init__(self):

        super().__init__()

        # Параметры волновода
        typeVolni = root.type_volni_entry.get()
        m = get(root.moda_m_entry)
        n = get(root.moda_n_entry)
        a = get(root.size_a_entry)  # Размер волновода по x
        b = get(root.size_b_entry)  # Размер волновода по y
        frequencyGGz = get(root.frequency_entry)

        self.type_polya = ""
        self.lenght_of_vector = 1
        self.XX = []
        self.YY = []
        self.ZZ = []

        cc = 3 * (10 ** 11)  # Скорость света 3 мм/с
        lambda_kr = (2 * a * b) / sqrt(pow(m * b, 2) + pow(n * a, 2))
        lambda1 = cc / (frequencyGGz * 1e9)

        # c = lambda1 * (m + n)  # Размер волновода по z
        c = get(root.size_c_entry)  # Размер волновода по z

        frequency = (frequencyGGz * 1e9)
        w = 2 * pi * frequency
        beta = (w / cc) * (sqrt(1 - pow((lambda1 / lambda_kr), 2)))
        kappa = sqrt((m * pi / a) ** 2 + (n * pi / b) ** 2)

        # Количество линий уровня
        # k = get(root.lines_entry)
        # h = a / k  # !!!!!! Шаг взять у пользователя

        f_kr = cc / lambda1
        T = 1 / f_kr

        Hm = 20
        Emax = 20
        Eps = 1
        mu = 1

        def selected(choice):
            if choice == "Не фиксировать" and self.visible:
                self.which_p_to_fix_value_label.grid_remove()
                self.which_p_to_fix_value.grid_remove()
                self.visible = False
            elif choice != "Не фиксировать" and not self.visible:
                self.which_p_to_fix_value_label.grid(row=4, column=0, sticky='w', padx=15)
                self.which_p_to_fix_value['validatecommand'] = (frameOfGraph.register(testVal), '%P', '%d')
                self.which_p_to_fix_value.grid(row=4, column=0, sticky='e', padx=15)
                self.visible = True

        def play():
            self.check = True
            while self.check:
                self.time += T / 40
                ReDrawQ()

        def time_minus():
            self.time -= T / 40
            ReDrawQ()

        def time_plus():
            self.time += T / 40
            ReDrawQ()

        def startPlay():
            if not self.check:
                thread = Thread(target=play, daemon=True)
                thread.start()

        def stopPlay():
            self.check = False

        def Formules(x, y, z, typeVolni, typePolya):
            if typeVolni == "TE":
                if typePolya == "E":
                    Xcomp = -(w * mu / (kappa ** 2)) * (n * pi / b) * Hm * cos(m * pi * x / a) * sin(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Ycomp = ((w * mu) / (kappa ** 2)) * (m * pi / a) * Hm * sin(m * pi * x / a) * cos(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Zcomp = x * y * z * 0
                    return Xcomp, Ycomp, Zcomp
                elif typePolya == "H":
                    Xcomp = -(beta / (kappa ** 2)) * (m * pi / a) * Hm * sin(m * pi * x / a) * cos(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Ycomp = -(beta / (kappa ** 2)) * (n * pi / b) * Hm * cos(m * pi * x / a) * sin(
                        n * pi * y / b) * sin(w * self.time - beta * z)
                    Zcomp = Hm * cos(m * pi * x / a) * cos(n * pi * y / b) * cos(w * self.time - beta * z)
                    return Xcomp, Ycomp, Zcomp
            else:
                if typePolya == "E":
                    Xcomp = (-(beta * m * pi) / ((kappa ** 2) * a)) * Emax * cos((m * pi * x) / a) * sin(
                        (n * pi * y) / b) * cos(w * self.time - beta * z)
                    Ycomp = (-(beta * n * pi) / ((kappa ** 2) * b)) * Emax * sin((m * pi * x) / a) * cos(
                        (n * pi * y) / b) * cos(w * self.time - beta * z)
                    Zcomp = Emax * sin((m * pi * x) / a) * sin((n * pi * y) / b) * sin(w * self.time - beta * z)
                    return Xcomp, Ycomp, Zcomp
                elif typePolya == "H":
                    Xcomp = -((w * Eps * n * pi) / ((kappa ** 2) * b)) * Emax * sin((m * pi * x) / a) * cos(
                        (n * pi * y) / b) * cos(w * self.time - beta * z)
                    Ycomp = ((w * Eps * m * pi) / ((kappa ** 2) * a)) * Emax * cos((m * pi * x) / a) * sin(
                        (n * pi * y) / b) * cos(w * self.time - beta * z)
                    Zcomp = 0 * x * y * z
                    return Xcomp, Ycomp, Zcomp

        def GetVolnovod(a, b, c):
            pogr = 1
            xx = [0 - pogr, a + pogr]
            yy = [0 - pogr, b + pogr]
            zz = [0 - pogr, c + pogr]

            X = np.array(
                [xx[0], xx[1], xx[1], xx[0], xx[0], xx[0], xx[1], xx[1], xx[1], xx[1], xx[1], xx[1], xx[0], xx[0],
                 xx[0], xx[0]])
            Y = np.array(
                [yy[0], yy[0], yy[0], yy[0], yy[0], yy[1], yy[1], yy[0], yy[1], yy[1], yy[0], yy[1], yy[1], yy[0],
                 yy[1], yy[1]])
            Z = np.array(
                [zz[0], zz[0], zz[1], zz[1], zz[0], zz[0], zz[0], zz[0], zz[0], zz[1], zz[1], zz[1], zz[1], zz[1],
                 zz[1], zz[0]])
            return X, Y, Z

        def getData(x, y, z, shag):
            fix_value = get(self.which_p_to_fix_value)
            if which_p_to_fix_entry.get() == "Не фиксировать":
                A, B, C = np.meshgrid(np.linspace(0, x, int(x / shag)),
                                      np.linspace(0, y, int(y / shag)),
                                      np.linspace(0, z, int(z / shag)))
                return A, B, C
            elif fix_value is None:
                return
            elif which_p_to_fix_entry.get() == "X":
                if fix_value > x or fix_value < 0:
                    tkinter.messagebox.showerror(title=None, message="Величина выходит за границы волновода!")
                    return
                A, B, C = np.meshgrid(fix_value,
                                      np.linspace(0, y, int(y / shag)),
                                      np.linspace(0, z, int(z / shag)))
                return A, B, C
            elif which_p_to_fix_entry.get() == "Y":
                if fix_value > y or fix_value < 0:
                    tkinter.messagebox.showerror(title=None, message="Величина выходит за границы волновода!")
                    return
                A, B, C = np.meshgrid(np.linspace(0, x, int(x / shag)),
                                      fix_value,
                                      np.linspace(0, z, int(z / shag)))
                return A, B, C
            elif which_p_to_fix_entry.get() == "Z":
                if fix_value > z or fix_value < 0:
                    tkinter.messagebox.showerror(title=None, message="Величина выходит за границы волновода!")
                    return
                A, B, C = np.meshgrid(np.linspace(0, x, int(x / shag)),
                                      np.linspace(0, y, int(y / shag)),
                                      fix_value)
                return A, B, C
            else:
                tkinter.messagebox.showerror(title=None, message="Неправильно введены данные!")
                return False

        def calculation():
            k = get(k_entry)
            if get(k_entry) == 0:
                tkinter.messagebox.showerror(title=None, message="Заполните поле количества векторов!")
                return
            if get(len_of_vector_entry) is None:
                tkinter.messagebox.showerror(title=None, message="Заполните поле длины векторов!")
                return
            if getData(a, b, c, min(a, b, c) / k) is not None:
                self.XX, self.YY, self.ZZ = getData(a, b, c, min(a, b, c) / k)

            self.type_polya = which_p_entry.get()
            self.lenght_of_vector = get(len_of_vector_entry)

            U, V, W = Formules(self.XX, self.YY, self.ZZ, typeVolni, self.type_polya)
            UVWmodule = (sqrt(U ** 2 + V ** 2 + W ** 2))

            if which_p_to_fix_entry.get() != "Не фиксировать":
                mappable = cm.ScalarMappable(cmap='turbo')
                self.colorBar = fig.colorbar(mappable, cax=cax, orientation='vertical')
            else:
                mappable = cm.ScalarMappable(cmap='turbo')
                mappable.set_array(np.linspace(0, np.max(UVWmodule), 10))
                self.colorBar = fig.colorbar(mappable, cax=cax, orientation='vertical')

            ReDrawQ()

        class FieldLine:
            """A Field Line."""

            def __init__(self, x):
                "Initializes the field line points 'x'."""
                self.x = x

            def plot(self, linewidth=None, linestyle='-',
                     startarrows=True, endarrows=True):
                """Plots the field line and arrows."""

                if linewidth is None:
                    linewidth = matplotlib.rcParams['lines.linewidth']

                x, y = zip(*self.x)
                pyplot.plot(x, y, '-k', linewidth=linewidth, linestyle=linestyle)

                n = int(len(x) / 2) if len(x) < 225 else 75
                if startarrows:
                    pyplot.arrow(x[n], y[n], (x[n + 1] - x[n]) / 100., (y[n + 1] - y[n]) / 100.,
                                 fc="k", ec="k",
                                 head_width=0.1 * linewidth, head_length=0.1 * linewidth)

                if len(x) < 225 or not endarrows:
                    return

                pyplot.arrow(x[-n], y[-n],
                             (x[-n + 1] - x[-n]) / 100., (y[-n + 1] - y[-n]) / 100.,
                             fc="k", ec="k",
                             head_width=0.1 * linewidth, head_length=0.1 * linewidth)

        def ReDrawQ():
            self.q.remove()
            const = 1e-10

            x, y, z = self.XX, self.YY, self.ZZ

            U, V, W = Formules(x, y, z, typeVolni, self.type_polya)

            uvw = np.vstack((U[np.newaxis], V[np.newaxis], W[np.newaxis]))
            for i in range(0, 3):
                if np.linalg.norm(uvw[i]) > 0:
                    uvw[i] += 1e-100
            U, V, W = uvw

            maximum = np.max(sqrt(U ** 2 + V ** 2 + W ** 2))

            U = U / maximum
            V = V / maximum
            W = W / maximum

            c = sqrt(U ** 2 + V ** 2 + W ** 2)
            c = (c.ravel() - c.min()) / c.ptp()
            c = np.concatenate((c, np.repeat(c, 2)))
            c = plt.cm.turbo(c)

            self.q = ax.quiver(x, z, y, U, W, V, cmap='turbo', colors=c, length=self.lenght_of_vector, normalize=False)
            #ax.scatter(U, W, V)

            canvas.draw()

        # Конфигурация окна
        self.title("График")
        self.w, self.h = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry("%dx%d" % (self.w, self.h))

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        frameOfGraph = customtkinter.CTkFrame(self, fg_color="gray", corner_radius=0)
        frameOfGraph.grid(row=0, column=0, sticky='nsew')

        frameOfGraph.grid_rowconfigure((0, 1, 2, 3, 4, 5), weight=1)
        frameOfGraph.grid_columnconfigure(0, weight=1)
        frameOfGraph.grid_columnconfigure(1, weight=2)

        ######1Column
        which_p = customtkinter.CTkLabel(frameOfGraph, text="Выберите какое поле построить: ", bg_color="transparent")
        which_p.grid(row=0, column=0, sticky='w', padx=15)

        which_p_entry = customtkinter.CTkOptionMenu(frameOfGraph, values=["E", "H"], fg_color="#565B5E",
                                                    button_color="#565B5E", button_hover_color="#464A4D",
                                                    dropdown_fg_color="#464A4D", dropdown_hover_color="#464A4D")
        which_p_entry.grid(row=0, column=0, sticky='e', padx=15)

        k_label = customtkinter \
            .CTkLabel(frameOfGraph, text="Введите количество векторов на оси\n(Влияет на скорость построения!) ",
                      bg_color="transparent", text_color="red") \
            .grid(row=1, column=0, sticky='w', padx=15)
        k_entry = customtkinter.CTkEntry(frameOfGraph, width=50)
        k_entry['validatecommand'] = (k_entry.register(testVal), '%P', '%d')
        k_entry.grid(row=1, column=0, sticky='e', padx=15)

        len_of_vector_label = customtkinter \
            .CTkLabel(frameOfGraph, text="Введите длину векторов! ",
                      bg_color="transparent") \
            .grid(row=2, column=0, sticky='w', padx=15)
        len_of_vector_entry = customtkinter.CTkEntry(frameOfGraph, width=50)
        len_of_vector_entry['validatecommand'] = (k_entry.register(testVal), '%P', '%d')
        len_of_vector_entry.grid(row=2, column=0, sticky='e', padx=15)

        which_p_to_fix = customtkinter.CTkLabel(frameOfGraph, text="Выберите какую ось зафиксировать: ",
                                                bg_color="transparent")
        which_p_to_fix.grid(row=3, column=0, sticky='w', padx=15)

        which_p_to_fix_entry = customtkinter.CTkOptionMenu(frameOfGraph, values=["Не фиксировать", "X", "Y", "Z"],
                                                           fg_color="#565B5E",
                                                           button_color="#565B5E", button_hover_color="#464A4D",
                                                           dropdown_fg_color="#464A4D", dropdown_hover_color="#464A4D",
                                                           command=selected)
        which_p_to_fix_entry.grid(row=3, column=0, sticky='e', padx=15)

        self.which_p_to_fix_value_label = customtkinter.CTkLabel(frameOfGraph,
                                                                 text="Значение, где зафиксировать ось: ",
                                                                 bg_color="transparent")
        self.which_p_to_fix_value = customtkinter.CTkEntry(frameOfGraph, width=50)

        frameOfChanging = customtkinter.CTkFrame(frameOfGraph, fg_color="gray", corner_radius=0)
        frameOfChanging.grid(row=5, column=0, rowspan=1, sticky='nsew')
        frameOfChanging.grid_rowconfigure((0, 1, 2), weight=1)
        frameOfChanging.grid_columnconfigure((0, 1), weight=1)

        CTkButton(frameOfChanging, text="Рассчитать", fg_color="#565B5E", hover_color="#464A4D", command=calculation) \
            .grid(row=0, columnspan=2)

        CTkButton(frameOfChanging, text="-", fg_color="#565B5E", hover_color="#464A4D", command=time_minus).grid(row=1,
                                                                                                                 column=0)
        CTkButton(frameOfChanging, text='+', fg_color="#565B5E", hover_color="#464A4D", command=time_plus).grid(row=1,
                                                                                                                column=1)
        CTkButton(frameOfChanging, text="Play", fg_color="#565B5E", hover_color="#464A4D", command=startPlay).grid(
            row=2, column=0)
        CTkButton(frameOfChanging, text="Stop", fg_color="#565B5E", hover_color="#464A4D", command=stopPlay).grid(row=2,
                                                                                                                  column=1)

        ######2Column

        fig = figure(facecolor="white")
        ax = fig.add_subplot(projection='3d')
        cax = fig.add_subplot([0.03, 0.5, 0.02, 0.47])
        lim = (0, max(a, b, c))
        plt.setp(ax, xlim=lim, ylim=lim, zlim=lim, xlabel="\nX", ylabel="\nZ", zlabel="\nY")
        X, Y, Z = GetVolnovod(a, b, c)
        ax.plot(X, Z, Y, color='r')

        self.q = ax.quiver(0, 0, 0, 0, 0, 0)

        self.colorBar = fig.colorbar(self.q, cax=cax, cmap="turbo")

        canvas = FigureCanvasTkAgg(fig, master=frameOfGraph)
        canvas.get_tk_widget().grid(row=0, rowspan=5, column=1, sticky="nsew")

        toolbarFrame = CTkFrame(master=frameOfGraph, corner_radius=0, fg_color='white')
        toolbarFrame.grid(row=5, column=1, sticky='nsew')
        toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)

        self.protocol("WM_DELETE_WINDOW", self.quit_clicked_close_plot)
        self.mainloop()

    def quit_clicked(self):
        if messagebox.askokcancel("Закрыть", "Действительно хотите закрыть окно?"):
            self.destroy()

    def quit_clicked_close_plot(self):
        self.check = False
        thread = Thread(target=self.quit_clicked, daemon=True)
        thread.start()


def quit_clicked():
    root.destroy()
    sys.exit()


root.protocol("WM_DELETE_WINDOW", quit_clicked)
root.mainloop()