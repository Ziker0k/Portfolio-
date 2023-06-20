import tkinter.messagebox

import customtkinter
from numpy import (sin, cos, tan, pi, exp, sqrt, arange, meshgrid, linspace)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.pyplot import (clf, gcf, subplots_adjust, figure)
from sys import exit
import os

customtkinter.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("dark-blue")  # Themes: "blue" (standard), "green", "dark-blue"

time = 0
class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        def testVal(inStr, acttyp):
            if acttyp == '1':  # insert
                if not inStr.isdigit():
                    return False
            return True

        def get(entry):
            value = entry.get()
            try:
                return int(value)
            except ValueError:
                return None

        def time_minus():
            global time
            if time > 0:
                time -= 1
                self.time_entry.configure(text=time)
                plotting()

        def time_plus():
            global time
            time += 1
            self.time_entry.configure(text=time)
            plotting()

        def plotting():
            global size_lines_entry
            global frequency
            global moda_n
            global moda_m
            global label
            global f_kr
            global time
            clf()

            if (len(self.moda_n_entry.get()) == 0) or (len(self.moda_m_entry.get()) == 0) or (len(self.size_a_entry.get()) == 0) or (len(self.size_b_entry.get()) == 0)\
                    or (len(self.frequency_entry.get()) == 0) or (len(self.lines_entry.get()) == 0):
                clf()
                tkinter.messagebox.showerror(title=None, message="Заполните все поля!")
                gcf().canvas.draw()
                return

            cc = 3e10  # Скорость света 6 см/с

            ff = get(self.frequency_entry)  # Частота 6 ГГц
            f = ff * 1e9  # Перевод в Гц

            w = 2 * pi * f
            lyam = cc / f
            hh = w / cc  # Волновое число в волноводе
            # Параметры волновода
            n = get(self.moda_n_entry)
            m = get(self.moda_m_entry)
            a = get(self.size_a_entry)  # Размер волновода по x
            b = get(self.size_b_entry)  # Размер волновода по y
            c = lyam  # Размер волновода по z

            # Шаг сетки
            h = 0.01
            # Время сечения
            tt = time  # Время
            t = tt / 1e12
            # Количество линий уровня
            k = get(self.lines_entry)

            kappa = sqrt(((pi * n / a) ** 2) + ((pi * m / b) ** 2))
            kappaX = pi * n / a
            kappaY = pi * n / b
            f_kr = (cc * kappa) / (2 * pi)

            # Алгоритм определения констант для проекций (где будет сделан срез)

            C1 = 0
            C2 = 0
            if n != 0 and m == 0:
                if (n % 2) != 0:
                    C1 = a / 2
                else:
                    C1 = a / (2 * n)
            elif n == 0 and m != 0:
                if (m % 2) != 0:
                    C2 = b / 2
                else:
                    C2 = b / (2 * m)

            ########################### TE Волны
            #### H
            def TE_H_XY(x, y):
                return ((abs(sin(kappaX * x))) / (abs(sin(kappaY * y)))) * cos(w * t + pi / 2)

            def TE_H_YZ(y, z):
                return (abs(sin(kappaY * y))) * cos(w * t - hh * z + pi / 2)

            def TE_H_XZ(x, z):
                return (abs(sin(kappaX * x))) * cos(w * t - hh * z + pi / 2)

            # E
            def TE_E_XY(x, y):
                return abs(cos(kappaY * y)) * abs(cos(kappaX * x)) * cos(w * t)

            # H10
            def TE10_H_XY(x, y):
                return abs(sin(kappaX * x)) * exp(hh * tan(w * t) * y)

            def TE10_H_XZ(x, z):
                return abs(sin(kappaX * x)) * cos(w * t - hh * z)

            def TE10_H_YZ(y, z):
                return exp(kappaX * (cos(kappaX * C1) / sin(kappaX * C1)) * y) * cos(w * t - hh * z)

            # E10
            1

            #################### H01 ######################
            def TE01_H_XY(x, y):
                return abs(sin(kappaY * y)) * exp(hh * tan(w * t) * x)

            def TE01_H_XZ(x, z):
                return exp(kappaY * (cos(kappaY * C2) / sin(kappaY * C2)) * x) * cos(w * t - hh * z)

            def TE01_H_YZ(y, z):
                return abs(sin(kappaY * y)) * cos(w * t - hh * z)

            #################### E01 #######################
            def TE01_E_XY(x, y):
                return (w / (kappaY * cc)) * abs(sin(kappaY * y)) * sin(w * t)

            ################################ TM ###################################
            ### H ###
            def TM_H_XY(x, y):
                return abs(sin(kappaX * x)) * abs(sin(kappaY * y)) * cos(w * t + pi / 2)

            ######### E #########
            def TM_E_XY(x, y):
                return (abs(cos(kappaY * y))) / (abs(cos(kappaX * x))) * cos(w * t)

            def TM_E_XZ(x, z):
                return (abs(cos(kappaX * x))) * cos(w * t - hh * z)

            def TM_E_YZ(y, z):
                return (abs(cos(kappaY * y))) * cos(w * t - hh * z)

            def makeData(b1, b2):
                a1 = arange(0, b1, h)
                a2 = arange(0, b2, h)
                a1grid, a2grid = meshgrid(a1, a2)
                return a1grid, a2grid

            if f > f_kr:
                XY = fig.add_subplot(3, 1, 1)
                YZ = fig.add_subplot(3, 1, 2)
                XZ = fig.add_subplot(3, 1, 3)

                fig.set_figwidth(5)
                fig.set_figheight(5)

                x1, y1 = makeData(a, b)
                y2, z2 = makeData(b, c)
                x3, z3 = makeData(a, c)
                if self.type_volni_entry.get() == "TE":
                    ###################### TE #########################
                    ########### XY ############
                    ##### H #####
                    if n != 0 and m != 0:
                        XY.contour(x1, y1, TE_H_XY(x1, y1), linspace(-1, 1, k), colors='b')
                    elif n != 0 and m == 0:
                        XY.contour(x1, y1, TE10_H_XY(x1, y1), linspace(-1, 1, k), colors='b')
                    elif n == 0 and m != 0:
                        XY.contour(x1, y1, TE01_H_XY(x1, y1), linspace(-1, 1, k), colors='b')
                    XY.set_xlabel('x')
                    XY.set_ylabel('y')

                    ####### E #########
                    if n != 0 and m != 0:
                        XY.contour(x1, y1, TE_E_XY(x1, y1), linspace(-1, 1, k), colors='r')
                    elif n != 0 and m == 0:
                        XY.contour(x1, y1, TE10_E_XY(x1, y1), linspace(-1, 1, k), colors='r')
                    elif n == 0 and m != 0:
                        XY.contour(x1, y1, TE01_E_XY(x1, y1), linspace(-1, 1, k), colors='r')
                    XY.set_xlabel('x')
                    XY.set_ylabel('y')
                    ############# YZ ##################
                    ###### H #######
                    if n != 0 and m != 0:
                        YZ.contour(z2, y2, TE_H_YZ(y2, z2), linspace(-1, 1, k), colors='b')
                    elif n != 0 and m == 0:
                        YZ.contour(z2, y2, TE10_H_YZ(y2, z2), linspace(-1, 1, k), colors='b')
                    elif n == 0 and m != 0:
                        YZ.contour(z2, y2, TE01_H_YZ(y2, z2), linspace(-1, 1, k), colors='b')
                    YZ.set_xlabel('z')
                    YZ.set_ylabel('y')

                    ############ XZ ###################

                    ########## H ###############
                    if n != 0 and m != 0:
                        XZ.contour(z3, x3, TE_H_XZ(x3, z3), linspace(-1, 1, k), colors='b')
                    elif n != 0 and m == 0:
                        XZ.contour(z3, x3, TE10_H_XZ(x3, z3), linspace(-1, 1, k), colors='b')
                    elif n == 0 and m != 0:
                        XZ.contour(z3, x3, TE01_H_XZ(x3, z3), linspace(-1, 1, k), colors='b')
                    XZ.set_xlabel('z')
                    XZ.set_ylabel('x')

                else:

                    ################################# TM ##################################
                    ################## XY #######################
                    ######### H ##########
                    if n != 0 and m != 0:
                        XY.contour(x1, y1, TM_H_XY(x1, y1), linspace(-1, 1, k), colors='b')
                        XY.set_xlabel('x')
                        XY.set_ylabel('y')
                        ######### E ###########
                        XY.contour(x1, y1, TM_E_XY(x1, y1), linspace(-1, 1, k), colors='r')
                        XY.set_xlabel('x')
                        XY.set_ylabel('y')

                        ################## YZ #######################
                        ######### E ###########
                        YZ.contour(z2, y2, TM_E_YZ(y2, z2), linspace(-1, 1, k), colors='r')
                        YZ.set_xlabel('z')
                        YZ.set_ylabel('y')

                        ################## XZ #######################
                        ########## E ##########
                        XZ.contour(z3, x3, TM_E_XZ(x3, z3), linspace(-1, 1, k), colors='r')
                        XZ.set_xlabel('z')
                        XZ.set_ylabel('x')
            else:
                tkinter.messagebox.showerror(title=None, message="Частота ниже критической!")
                self.label.configure(text="Критическая частота(ГГц): " + str(f_kr / 1e9))
                return
            subplots_adjust(wspace=0.5, hspace=0.5)
            gcf().canvas.draw()
            self.label.configure(text="Критическая частота(ГГц): " + str(f_kr / 1e9))

        # configure window
        self.title("Волновод")
        self.geometry(f"{1100}x{600}")

        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0, 1, 2), weight=1)

        # create sidebar frame with widgets
        self.sidebar_frame = customtkinter.CTkFrame(self, width=300, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(7, weight=1)

        self.logo_label = customtkinter.CTkLabel(self.sidebar_frame, text="Характеристики волновода", font=customtkinter.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, columnspan=2, padx=20, pady=(20, 10))

        #row1
        self.label_n = customtkinter.CTkLabel(self.sidebar_frame, text="Выберите тип волны:")
        self.label_n.grid(row=1, column=0, sticky="w", padx=15)
        self.type_volni_entry = customtkinter.CTkOptionMenu(self.sidebar_frame, values=["TE", "TM"], fg_color="#565B5E", button_color="#565B5E", button_hover_color="#464A4D", dropdown_fg_color="#464A4D", dropdown_hover_color="#464A4D")
        self.type_volni_entry.grid(row=1, column=1, padx=15)

        #row2
        self.moda_label = customtkinter.CTkLabel(self.sidebar_frame, text="Выберите моду:")
        self.moda_label.grid(row=2, column=0, sticky="w", padx=15, pady=(15, 0))

        #row3
        self.moda_n = customtkinter.CTkLabel(self.sidebar_frame, text="n: ")
        self.moda_n.grid(row=3, column=0, sticky="w", padx=15)
        self.moda_n_entry = customtkinter.CTkEntry(self.sidebar_frame, width=50)
        self.moda_n_entry['validatecommand'] = (self.moda_n_entry.register(testVal), '%P', '%d')
        self.moda_n_entry.grid(row=3, column=0, padx=(30, 0), sticky="w")
        self.moda_m = customtkinter.CTkLabel(self.sidebar_frame, text="m: ")
        self.moda_m.grid(row=3, column=1, sticky="w", padx=15)
        self.moda_m_entry = customtkinter.CTkEntry(self.sidebar_frame, width=50)
        self.moda_m_entry['validatecommand'] = (self.moda_m_entry.register(testVal), '%P', '%d')
        self.moda_m_entry.grid(row=3, column=1, padx=(30, 0), sticky="w")

        #row4
        self.frequency_label = customtkinter.CTkLabel(self.sidebar_frame, text="Введите частоту(ГГц): ")
        self.frequency_label.grid(row=4, column=0, padx=15, pady=15, sticky="w")
        self.frequency_entry = customtkinter.CTkEntry(self.sidebar_frame)
        self.frequency_entry['validatecommand'] = (self.frequency_entry.register(testVal), '%P', '%d')
        self.frequency_entry.grid(row=4, column=1, sticky="w", padx=15)

        self.label = customtkinter.CTkLabel(self.sidebar_frame, text="")
        self.label.grid(row=5, columnspan=2, padx=15, pady=15, sticky="SW")

        #row5
        self.size_label = customtkinter.CTkLabel(self.sidebar_frame, text="Введите размеры волновода(см):")
        self.size_label.grid(row=6, column=0, sticky="w", padx=15)
        self.size_a_label = customtkinter.CTkLabel(self.sidebar_frame, text="a:")
        self.size_a_label.grid(row=6, column=1, sticky="w", padx=15)
        self.size_a_entry = customtkinter.CTkEntry(self.sidebar_frame, width=50)
        self.size_a_entry['validatecommand'] = (self.size_a_entry.register(testVal), '%P', '%d')
        self.size_a_entry.grid(row=6, column=1, padx=(30, 0), sticky="w")

        #row6
        self.size_b_label = customtkinter.CTkLabel(self.sidebar_frame, text="b:")
        self.size_b_label.grid(row=7, column=1, sticky="w", padx=15)
        self.size_b_entry = customtkinter.CTkEntry(self.sidebar_frame, width=50)
        self.size_b_entry['validatecommand'] = (self.size_b_entry.register(testVal), '%P', '%d')
        self.size_b_entry.grid(row=7, column=1, padx=(30, 0), sticky="w")

        #row7
        self.lines_label = customtkinter.CTkLabel(self.sidebar_frame, text="Введите количество линий уровня:", anchor="w")
        self.lines_label.grid(row=8, column=0, padx=15)
        self.lines_entry = customtkinter.CTkEntry(self.sidebar_frame, width=50)
        self.lines_entry['validatecommand'] = (self.lines_entry.register(testVal), '%P', '%d')
        self.lines_entry.grid(row=8, column=1, sticky="w", padx=15)

        #row8
        self.time_label = customtkinter.CTkLabel(self.sidebar_frame, text="Время:")
        self.time_label.grid(row=9, column=0, sticky="w", padx=15, pady=15)

        self.time_minus = customtkinter.CTkButton(self.sidebar_frame, text="-", command=time_minus, width=50, fg_color="#565B5E", hover_color="#464A4D")
        self.time_minus.grid(row=9, column=1, sticky="w", padx=15)
        self.time_entry = customtkinter.CTkLabel(self.sidebar_frame, text=str(0))
        self.time_entry.grid(row=9, column=1)
        self.time_plus = customtkinter.CTkButton(self.sidebar_frame, text="+", command=time_plus, width=50, fg_color="#565B5E", hover_color="#464A4D")
        self.time_plus.grid(row=9, column=1, sticky="e", padx=15)

        #row10
        self.plotting_button = customtkinter.CTkButton(self.sidebar_frame, text="Построить", command=plotting,
                                                       fg_color="#565B5E", hover_color="#464A4D")
        self.plotting_button.grid(row=10, columnspan=2, pady=15)

        #2ond frame
        self.frame2 = customtkinter.CTkFrame(self, fg_color="lightgray", corner_radius=0)
        self.minsize(1100, 600)
        self.frame2.grid(row=0, columnspan=3, column=1, rowspan=4, sticky="nsew")
        self.frame2.grid_rowconfigure(1, weight=1)

        label_color_blue = customtkinter.CTkLabel(self.frame2, text="Н-синий", bg_color="transparent", text_color="blue")
        label_color_blue.grid(row=1, column=0, columnspan=1)
        label_color_red = customtkinter.CTkLabel(self.frame2, text="E-красный", bg_color="transparent", text_color="red")
        label_color_red.grid(row=2, column=0, columnspan=1, pady=(0, 15))

        fig = figure(figsize=(5 , 5), dpi=100, facecolor="lightgray")
        canvas = FigureCanvasTkAgg(fig, master=self.frame2)
        canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')

        self.minsize(self.sidebar_frame.winfo_width()+600, 500)
        self.maxsize(self.winfo_width(), self.winfo_height())



        # set default values
        #self.slider_2.configure(command=self.progressbar_3.set)

    def change_appearance_mode_event(self, new_appearance_mode: str):
        customtkinter.set_appearance_mode(new_appearance_mode)

    def change_scaling_event(self, new_scaling: str):
        new_scaling_float = int(new_scaling.replace("%", "")) / 100
        customtkinter.set_widget_scaling(new_scaling_float)


if __name__ == "__main__":
    app = App()
    app.mainloop()