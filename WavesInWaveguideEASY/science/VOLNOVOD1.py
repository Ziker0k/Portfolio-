import sys
from threading import Thread

import customtkinter
import matplotlib
import tkinter.messagebox
import numpy as np
import time as tm
from numpy import (sin, cos, pi, linspace, arange, sqrt)

matplotlib.use('TkAgg')
from matplotlib import pyplot as plt, cm, gridspec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.pyplot import (clf, gcf, subplots_adjust, figure)
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


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


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()
        self.check = False
        self.shag = 0
        self.time = 1
        self.n = None
        self.m = None
        self.a = None  # Размер волновода по x
        self.b = None  # Размер волновода по y
        self.frequencyGGz = None
        self.type_wave = None
        self.cc = 3 * (10 ** 11)  # Скорость света 3 мм/с
        self.shag = None
        Eps = 1
        mu = 1

        def time_minus():
            if self.time > 0:
                self.time -= self.shag
                #self.time_entry.configure(text=self.time)
                plotting()

        def time_plus():
            self.time += self.shag
            #self.time_entry.configure(text=self.time)
            plotting()

        def play():
            self.check = True
            while self.check:
                self.time += self.shag
                plotting()

        def startPlay():
            if not self.check:
                thread = Thread(target=play, daemon=True)
                thread.start()

        def stopPlay():
            self.check = False

        def plotting():
            if (len(self.moda_n_entry.get()) == 0) or (len(self.moda_m_entry.get()) == 0) or (
                    len(self.size_a_entry.get()) == 0) or (len(self.size_b_entry.get()) == 0) \
                    or (len(self.frequency_entry.get()) == 0):
                tkinter.messagebox.showerror(title=None, message="Заполните все поля!")
                gcf().canvas.draw()
                return

            # Параметры волновода
            n = get(self.moda_n_entry)
            m = get(self.moda_m_entry)
            a = get(self.size_a_entry)  # Размер волновода по x
            b = get(self.size_b_entry)  # Размер волновода по y
            frequencyGGz = get(self.frequency_entry)
            type_wave = self.type_wave_entry.get()

            if type_wave == "TM" and (m == 0 or n == 0):
                tkinter.messagebox.showerror(title=None, message="При моде TM m, n > 0!")
                gcf().canvas.draw()
                return

            cc = 3 * (10 ** 11)  # Скорость света 3 мм/с
            lambda_kr = (2 * a * b) / np.sqrt(pow(m * b, 2) + pow(n * a, 2))
            lambda1 = cc / (frequencyGGz * 1e9)

            frequency = (frequencyGGz * 1e9)
            w = 2 * np.pi * frequency
            beta = (w / cc) * (np.sqrt(1 - pow((lambda1 / lambda_kr), 2)))
            kappa = np.sqrt((m * np.pi / a) ** 2 + (n * np.pi / b) ** 2)
            c = lambda1

            f_kr = cc / lambda_kr
            T = 1 / frequency

            self.shag = T / 40

            Hm = 20
            Emax = 10000
            Eps = 1
            mu = 1

            # Алгоритм определения констант для проекций (где будет сделан срез)

            C1 = a
            C2 = b
            if (m % 2) != 0:
                C1 = a / 2
            else:
                C1 = a / (2 * (m or 1))
            if (n % 2) != 0:
                C2 = b / 2
            else:
                C2 = b / (2 * (n or 1))

            print(C1, "C1")
            print(C2, "C2")

            def Formules(x, y, z, typeVolni, typePolya):
                if typeVolni == "TE":
                    if typePolya == "E":
                        Xcomp = -cos(m * pi * x / a) * sin(n * pi * y / b) * sin(w * self.time - beta * z)
                        Ycomp = sin(m * pi * x / a) * cos(n * pi * y / b) * sin(w * self.time - beta * z)
                        Zcomp = x * y * z * 0
                        return Xcomp, Ycomp, Zcomp
                    elif typePolya == "H":
                        Xcomp = -sin(m * pi * x / a) * cos(n * pi * y / b) * sin(w * self.time - beta * z)
                        Ycomp = -cos(m * pi * x / a) * sin(n * pi * y / b) * sin(w * self.time - beta * z)
                        Zcomp = cos(m * pi * x / a) * cos(n * pi * y / b) * cos(w * self.time - beta * z)
                        return Xcomp, Ycomp, Zcomp
                else:
                    if typePolya == "E":
                        Xcomp = cos((m * pi * x) / a) * sin((n * pi * y) / b) * sin(w * self.time - beta * z)
                        Ycomp = sin((m * pi * x) / a) * cos((n * pi * y) / b) * sin(w * self.time - beta * z)
                        Zcomp = sin((m * pi * x) / a) * sin((n * pi * y) / b) * cos(w * self.time - beta * z)
                        return Xcomp, Ycomp, Zcomp
                    elif typePolya == "H":
                        Xcomp = -sin((m * pi * x) / a) * cos((n * pi * y) / b) * sin(w * self.time - beta * z)
                        Ycomp = cos((m * pi * x) / a) * sin((n * pi * y) / b) * sin(w * self.time - beta * z)
                        Zcomp = 0 * x * y * z
                        return Xcomp, Ycomp, Zcomp

            if frequency > f_kr:

                self.XY.clear()
                self.YZ.clear()
                self.XZ.clear()

                x = np.linspace(0, a, int(a))
                y = np.linspace(0, b, int(b))
                z = np.linspace(0, c, int(c))

                X, Y, Z = np.meshgrid(x, y, z)

                X1, Y1 = np.meshgrid(x, y)
                Z2, Y2 = np.meshgrid(z, y)
                Z3, X3 = np.meshgrid(z, x)

                EXmax, EYmax, EZmax = Formules(X, Y, Z, type_wave, "E")
                HXmax, HYmax, HZmax = Formules(X, Y, Z, type_wave, "H")

                EX1, EY1, EZ1 = Formules(X1, Y1, 0, type_wave, "E")
                HX1, HY1, HZ1 = Formules(X1, Y1, 0, type_wave, "H")

                maxE = np.max((sqrt(EXmax ** 2 + EYmax ** 2 + EZmax ** 2)))
                maxH = np.max((sqrt(HXmax ** 2 + HYmax ** 2 + HZmax ** 2)))

                module = sqrt(EX1 ** 2 + EY1 ** 2 + EZ1 ** 2) / (maxE)
                module2 = sqrt(HX1 ** 2 + HY1 ** 2 + HZ1 ** 2) / (maxH)


                density = (0.2*(m+1), 0.2*(n+1))

                self.XY.set_xlabel("X")
                self.XY.set_ylabel("Y")
                self.XY.streamplot(X1, Y1, EX1, EY1, linewidth=(module * 2), density=density, color='r',
                                   broken_streamlines=False,
                                   arrowstyle="-|>")
                self.XY.streamplot(X1, Y1, HX1, HY1, linewidth=(module2 * 2), density=density, color='b',
                                   broken_streamlines=False, arrowstyle='->')

                EX2, EY2, EZ2 = Formules(C1, Y2, Z2, type_wave, "E")
                HX2, HY2, HZ2 = Formules(C1, Y2, Z2, type_wave, "H")

                module3 = sqrt(EX2 ** 2 + EY2 ** 2 + EZ2 ** 2) / (maxE)
                module4 = sqrt(HX2 ** 2 + HY2 ** 2 + HZ2 ** 2) / (maxH)
                #
                self.YZ.set_xlabel("Z")
                self.YZ.set_ylabel("Y")
                self.YZ.streamplot(Z2, Y2, EZ2, EY2, density=density, color='r',
                                   broken_streamlines=False,
                                   linewidth=(module3*2),
                                   arrowstyle="-|>")
                self.YZ.streamplot(Z2, Y2, HZ2, HY2, density=density, color='b', linewidth=module4*2,arrowstyle="->",
                                   broken_streamlines=False)

                EX3, EY3, EZ3 = Formules(X3, C2, Z3, type_wave, "E")
                HX3, HY3, HZ3 = Formules(X3, C2, Z3, type_wave, "H")

                module5 = sqrt(EX3 ** 2 + EY3 ** 2 + EZ3 ** 2) / (maxE)
                module6 = sqrt(HX3 ** 2 + HY3 ** 2 + HZ3 ** 2) / (maxH)

                self.XZ.set_xlabel("Z")
                self.XZ.set_ylabel("X")
                self.XZ.streamplot(Z3, X3, EZ3, EX3, density=density, color='r',
                                   broken_streamlines=False, linewidth=module5*2,
                                   arrowstyle="-|>")
                self.XZ.streamplot(Z3, X3, HZ3, HX3, density=density, color='b', linewidth=module6*2,
                                   broken_streamlines=False, arrowstyle="->")
            else:
                tkinter.messagebox.showerror(title=None, message="Частота ниже критической!")
                self.label.configure(text="Критическая частота(ГГц): " + str(f_kr / 1e9))
                return
            self.canvas.draw()
            self.label.configure(text="Критическая частота(ГГц): " + str(f_kr / 1e9))

        # configure window
        self.title("Волновод")
        self.geometry(f"{1100}x{580}")

        # configure grid layout (2x1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=3)
        self.grid_rowconfigure(0, weight=1)

        # sidebar grid configure
        self.sidebar_frame = customtkinter.CTkFrame(self, width=300, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure((0, 1, 2, 3, 4, 5, 6, 7, 8, 9), weight=1)
        self.sidebar_frame.grid_columnconfigure((0, 1), weight=1)

        # row0
        self.logo_label = customtkinter.CTkLabel(self.sidebar_frame, text="Характеристики волновода",
                                                 font=customtkinter.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, columnspan=2, padx=20, pady=(20, 10))

        # row1
        self.type_wave_label = customtkinter.CTkLabel(self.sidebar_frame, text="Выберите тип волны:")
        self.type_wave_label.grid(row=1, column=0, sticky="w", padx=15)
        self.type_wave_entry = customtkinter.CTkOptionMenu(self.sidebar_frame, values=["TE", "TM"], fg_color="#565B5E",
                                                           button_color="#565B5E", button_hover_color="#464A4D",
                                                           dropdown_fg_color="#464A4D", dropdown_hover_color="#464A4D")
        self.type_wave_entry.grid(row=1, column=1, padx=15)

        # row2
        self.moda_label = customtkinter.CTkLabel(self.sidebar_frame, text="Выберите моду:")
        self.moda_label.grid(row=2, column=0, sticky="w", padx=15, pady=(15, 0))

        # row3
        self.moda_m = customtkinter.CTkLabel(self.sidebar_frame, text="m: ")
        self.moda_m.grid(row=3, column=0, sticky="w", padx=15)
        self.moda_m_entry = customtkinter.CTkEntry(self.sidebar_frame, width=50)
        self.moda_m_entry['validatecommand'] = (self.moda_m_entry.register(testVal), '%P', '%d')
        self.moda_m_entry.grid(row=3, column=0, padx=(30, 0), sticky="w")
        self.moda_n = customtkinter.CTkLabel(self.sidebar_frame, text="n: ")
        self.moda_n.grid(row=3, column=1, sticky="w", padx=15)
        self.moda_n_entry = customtkinter.CTkEntry(self.sidebar_frame, width=50)
        self.moda_n_entry['validatecommand'] = (self.moda_n_entry.register(testVal), '%P', '%d')
        self.moda_n_entry.grid(row=3, column=1, padx=(30, 0), sticky="w")

        # row4
        self.frequency_label = customtkinter.CTkLabel(self.sidebar_frame, text="Введите частоту(ГГц): ")
        self.frequency_label.grid(row=4, column=0, padx=15, pady=15, sticky="w")
        self.frequency_entry = customtkinter.CTkEntry(self.sidebar_frame)
        self.frequency_entry['validatecommand'] = (self.frequency_entry.register(testVal), '%P', '%d')
        self.frequency_entry.grid(row=4, column=1, sticky="w", padx=15)

        self.label = customtkinter.CTkLabel(self.sidebar_frame, text="")
        self.label.grid(row=5, columnspan=2, padx=15, pady=15, sticky="SW")

        # row6
        self.size_label = customtkinter.CTkLabel(self.sidebar_frame, text="Введите размеры волновода(мм):")
        self.size_label.grid(row=6, column=0, sticky="w", padx=15)

        # row7
        self.size_frame = customtkinter.CTkFrame(self.sidebar_frame, corner_radius=0)
        self.size_frame.grid(row=7, columnspan=2, sticky="nsew")
        self.size_frame.grid_columnconfigure((0, 1), weight=1)
        self.size_frame.grid_rowconfigure((0, 1), weight=1)

        self.size_a_label = customtkinter.CTkLabel(self.size_frame, text="a:")
        self.size_a_label.grid(row=0, column=0, sticky="nsew", padx=15)
        self.size_a_entry = customtkinter.CTkEntry(self.size_frame, width=50)
        self.size_a_entry['validatecommand'] = (self.size_a_entry.register(testVal), '%P', '%d')
        self.size_a_entry.grid(row=1, column=0, padx=15, pady=(0, 15), sticky="nsew")

        self.size_b_label = customtkinter.CTkLabel(self.size_frame, text="b:")
        self.size_b_label.grid(row=0, column=1, sticky="nsew", padx=15)
        self.size_b_entry = customtkinter.CTkEntry(self.size_frame, width=50)
        self.size_b_entry['validatecommand'] = (self.size_a_entry.register(testVal), '%P', '%d')
        self.size_b_entry.grid(row=1, column=1, padx=15, pady=(0, 15), sticky="nsew")

        frameOfChanging = customtkinter.CTkFrame(self.sidebar_frame, corner_radius=0)
        frameOfChanging.grid(row=8, column=0, columnspan=2, rowspan=1, sticky='nsew')
        frameOfChanging.grid_rowconfigure((0, 1, 2), weight=1)
        frameOfChanging.grid_columnconfigure((0, 1), weight=1)

        customtkinter.CTkButton(frameOfChanging, text="Рассчитать", fg_color="#565B5E", hover_color="#464A4D", command=plotting) \
            .grid(row=0, columnspan=2)

        customtkinter.CTkButton(frameOfChanging, text="-", fg_color="#565B5E", hover_color="#464A4D", command=time_minus).grid(row=1,
                                                                                                                 column=0)
        customtkinter.CTkButton(frameOfChanging, text='+', fg_color="#565B5E", hover_color="#464A4D", command=time_plus).grid(row=1,
                                                                                                                column=1)
        customtkinter.CTkButton(frameOfChanging, text="Play", fg_color="#565B5E", hover_color="#464A4D", command=startPlay).grid(
            row=2, column=0)
        customtkinter.CTkButton(frameOfChanging, text="Stop", fg_color="#565B5E", hover_color="#464A4D", command=stopPlay).grid(row=2,
                                                                                                                  column=1)

        # right frame
        self.right_frame = customtkinter.CTkFrame(self, fg_color="white", corner_radius=0)
        self.right_frame.grid(row=0, column=1, sticky='nsew')

        self.right_frame.grid_rowconfigure((0, 1, 2, 3), weight=2)
        self.right_frame.grid_rowconfigure((4, 5), weight=0)
        self.right_frame.grid_columnconfigure(0, weight=1)

        self.fig = figure(facecolor="white", constrained_layout=True, dpi=100)
        self.spec = gridspec.GridSpec(ncols=5, nrows=4, figure=self.fig)
        self.XY = self.fig.add_subplot(self.spec[:1, :2])
        self.XY.set_xlabel("X")
        self.XY.set_ylabel("Y")
        self.XY.set_xlim(0, 24)
        self.XY.set_ylim(0, 11)
        self.YZ = self.fig.add_subplot(self.spec[1:2, :])
        self.YZ.set_xlabel("Z")
        self.YZ.set_ylabel("Y")
        self.YZ.set_xlim(0, 50)
        self.YZ.set_ylim(0, 11)
        self.XZ = self.fig.add_subplot(self.spec[2:4, :])
        self.XZ.set_xlabel("Z")
        self.XZ.set_ylabel("X")
        self.XZ.set_xlim(0, 50)
        self.XZ.set_ylim(0, 24)
        self.XYplot = None
        self.YZplot = None
        self.ZYplot = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
        self.canvas.get_tk_widget().grid(row=0, rowspan=4, column=0, sticky="nsew", padx=30, pady=30)

        self.toolbarFrame = customtkinter.CTkFrame(master=self.right_frame, corner_radius=0, fg_color='white')
        self.toolbarFrame.grid(row=5, column=0, sticky='nsew')
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)

        label_color_blue = customtkinter.CTkLabel(self.right_frame, text="Н-синий", bg_color="transparent",
                                                  text_color="blue")
        label_color_blue.grid(row=4, column=0, sticky="e", padx=15)
        label_color_red = customtkinter.CTkLabel(self.right_frame, text="E-красный", bg_color="transparent",
                                                 text_color="red")
        label_color_red.grid(row=4, column=0, sticky='w', padx=15)

        self.protocol("WM_DELETE_WINDOW", self.quit_clicked)

    def quit_clicked(self):
        self.destroy()
        sys.exit()


if __name__ == "__main__":
    app = App()
    app.mainloop()
