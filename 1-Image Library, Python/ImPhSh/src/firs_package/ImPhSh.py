from PIL import Image,ImageDraw, ImageFilter
import numpy as np
import random
import math as m
from random import randint
import os



def TR(path):
    apo = Image.open(path)
    apoBW = apo.rotate(180)
    apoN = apo

    size1 = apo.size #Массив размера изображения
    rx = np.random.random() * (size1[0] / 2.)
    r2x = size1[0] - rx
    while(rx <= 1/3 * size1[0]):
        rx = np.random.random() * (size1[0] / 2.)

        r2x = size1[0] - rx



    k = size1[1]/size1[0]

    ry = k * rx

    r2y = size1[1] - ry


    #Фильтры
    draw = ImageDraw.Draw(apoBW)
    pix = apoBW.load()
    for i in range(size1[0]):
        for j in range(size1[1]):
            a = pix[i, j][0]
            b = pix[i, j][1]
            c = pix[i, j][2]
            S = (a + b + c) // 3
            draw.point((i, j), (S, S, S))
    del draw
    draw = ImageDraw.Draw(apoN)
    pix2 = apoN.load()
    for i in range(size1[0]):
        for j in range(size1[1]):
            a = pix2[i, j][0]
            b = pix2[i, j][1]
            c = pix2[i, j][2]
            draw.point((i, j), (255 - a, 255 - b, 255 - c))
    del draw        

    mask_im = Image.new("L", apo.size, 0)
    draw = ImageDraw.Draw(mask_im)
    draw.polygon(
        xy=(
            (0, 0),
            (rx, 0),
            (0, ry)
            ), fill=255
        )
    draw.polygon(
        xy=(
            (size1[0], size1[1]),
            (size1[0],r2y),
            (r2x, size1[1])
            ), fill=255
        )
    del draw
    #mask_im.save('circle.jpg', quality=95)


    apoN.paste(apoBW, (0, 0), mask_im)


    apoN.show()
    
    
def Photo_to_Sim(path, T):
    def rand(X, a): # Создаю массив рандомных точек
    
        while(X[0] <= 1. / 10. or X[0] >= 3. / 10.): #1/10 3/10, 
            X[0] = random.random()
        X[0] = X[0] * a.size[0]
    
        while(X[1] <= 7. / 10. or X[1] >= 9. / 10.):# 7/10 9/10
            X[1] = random.random()
        X[1] = X[1] * a.size[1]
    
        return X
    
    def sort(X): # Сортирую массив Х1 Х2 Х1 меньшее Х2 большее
        if (X[0] > X[1]):
            a = X[0]
            X[0] = X[1]
            X[1] = a
            return X
        if (X[0] < X[1]):
            return X
    
    l = T
    img = Image.open(path)
    img_draw = ImageDraw.Draw(img)

    X = list()

    X.append(random.random())
    X.append(random.random())

    X = sort(rand(X, img))
    N = 10.
    y = np.arange(0, img.size[1], 1)
    x1 = list()
    x2 = list()
    N = float(l)
    A = 50
    consty = (2* m.pi )/img.size[1]


    for i in y:
        x1.append(int(m.sin(i / N) * N + X[0])) # задаем значения для синуса
        x2.append((int(m.sin(i / N) * N + X[1])))

    img_draw.rectangle((0, 0, X[0], img.size[1] ), fill = (0, 0, 0, 0))
    img_draw.rectangle((X[1], 0, img.size[0], img.size[1] ), fill = (0, 0, 0, 0))
    

    pixdata = img.load()
    for yy in range(int(img.size[1])):
        if m.sin(yy*consty*N) > 0:
            for xx in range(int(X[0]), int(X[1] + A*m.sin(yy*consty*N)+1),1):
                pixdata[int(xx - A*m.sin(yy*consty*N)),yy] = pixdata[xx, yy]
        if m.sin(yy*consty*N) <= 0:
            for xx in range(int(X[1]), int(X[0] + A*m.sin(yy*consty*N)), -1):
                pixdata[int(xx - A*m.sin(yy*consty*N)),yy] = pixdata[xx, yy]
            
    img.show()
def Rand_Rectangle(path):
    im = Image.open(path)
    draw = ImageDraw.Draw(im)


    w=im.size[0] # запоминаем размеры изображения в пикселях
    h=im.size[1]

    pixelm=im.load()

    for i in range(w): # применяем фильтры
        for j in range(h):
            a=pixelm[i,j][0]
            b=pixelm[i,j][1]
            c=pixelm[i,j][2]
            S=(a+b+c)//3
            if (((i>=(w/4)) and (i<=(w/2 -1)))or((i>=((3*w)/4)) and(i<=(w-1)))):  # диапазон ширины 2 и 6 части  
                draw.point((i,j),(S,S,S))
              
            else:
                draw.point((i,j),(255-a, 255-b, 255-c))


    imn = im

    imn2 = imn.crop((w/4-1, 0, w/2-1, h/2-1))
    imn6 = imn.crop((3*w/4-1, 0, w-1, h/2-1))
    imn4 = imn.crop((w/4-1, h/2, w/2-1, h))
    imn8 = imn.crop((3*w/4-1, h/2, w-1, h))


    imn.paste(imn4,(480,0))
    imn.paste(imn8,(1440,0))
    imn.paste(imn2,(480,540))
    imn.paste(imn6,(1440,540))
    imn.show()


def Circule(path):
    img = Image.open(path)

    w, h = img.size

    wCenter  = w // 2
    hCenter = h // 2

    x = (w - h) // 2
    img_cropped = img.crop((x, 0, x + h, h))


    size =(h - 100, h - 100)
    mask = Image.new('L', size, 0)
    draw = ImageDraw.Draw(mask)

    draw.ellipse((0, 0) + size, fill = 255)

    img_cropped = img_cropped.resize(size)


    img_cropped1 = img_cropped
    img_cropped = img_cropped.convert('RGBA')
    pixdata = img_cropped.load()
    for yy in range(img_cropped.size[1]):
        for xx in range(img_cropped.size[0]):
            pixdata[xx,yy] = ( 0, 0 ,0)
    img_cropped.paste(img_cropped1, (0, 0), mask)

    img2=Image.new("RGBA",(w, h),( 0, 0 ,0))
    img2.paste(img_cropped,((wCenter - size[0] // 2, hCenter - size[0] // 2)))
    img2.show()
  
def Lines(path):
    ROW_COUNT = 12
    #----------------------------------------------------------------------------#
    
    #---------------------------Вспомогательные функции--------------------------#
    
    #Циклический сдвиг вправо
    def circular_shift(lst, shift):
        return lst[-shift:] + lst[:-shift]
    
    #Наибольший общий делитель
    def gcd(a,b):
        if b == 0:
            return a
        return gcd(b, a%b)
    
    #Негатив
    def get_negative(img):
        width = img.size[0]
        height = img.size[1]
        draw = ImageDraw.Draw(img)
        pix = img.load()
    
        #В цикле достаем значения r, g и b 
        #каждого пикселя и находим обратное им занченние
        for i in range(width):
            for j in range(height):
                r = pix[i, j][0]
                g = pix[i, j][1]
                b = pix[i, j][2]
                draw.point((i, j), (255 - r, 255 - g, 255 - b))
        return img
    
    #----------------------------------------------------------------------------#
    
    #-------------------------------Задание №1-----------------------------------#
    
    img = Image.open(path)
    

    
    width = img.size[0]
    height = img.size[1]
    
    #Делим изображение на 12 равных рядов
    row_height = int(height / ROW_COUNT)
    
    image_rows = []
    

    
    #Заносим ряды в список, чтобы случайно не перезаписать
    for i in range(ROW_COUNT):
        image_rows.append(
            img.crop((0, i*row_height, width, (i+1)*row_height))
        )
    
    #Меняем ряды местами и применям фильтры
    for i in range(ROW_COUNT):
        new_i = i
        if (i % 2 == 0):
            image_rows[i] = get_negative(image_rows[i])
            new_i = i + 2 if i + 2 < ROW_COUNT else 0
        else:
            image_rows[i] = image_rows[i].filter(ImageFilter.EDGE_ENHANCE_MORE)
            new_i = i - 2 if i - 2 >= 0 else ROW_COUNT - 1
        img.paste(
            image_rows[i],
            (0, new_i*row_height, width, (new_i+1)*row_height)
        )
    
    #Выдаем результат
    img.show()

def Carousel(path, border_scale, step_count):
    ROW_COUNT = 12
    def circular_shift(lst, shift):
        return lst[-shift:] + lst[:-shift]
    
    def gcd(a,b):
        if b == 0:
            return a
        return gcd(b, a%b)

    
    
    def get_new_image(orig_img, width, height, border_width, step_count):
    
        img = orig_img.copy()
    
        #Коробки для рамок
        border_boxes = [
            (border_width, height - border_width, width, height),           #Низ
            (width - border_width, 0, width, height - border_width),        #Право
            (0, 0, width - border_width, border_width),                     #Верх
            (0, border_width, border_width, height)                         #Лево
        ]
    
        #Вырезаем каждую сторону рамки и приводим их к одному виду
        borders = []
        for i in range(len(border_boxes)):
            borders.append(img.crop(border_boxes[i]))
            for j in range(i):
                borders[i] = borders[i].transpose(Image.ROTATE_270)
    
        #Подсчитываем НОД, чтобы получить прямоугольники одинаковых размеров
        gcd_result = gcd(borders[0].size[0], borders[1].size[0])
        #print("НОД: ", gcd_result)
    
        #Длина списка borders
        borders_length = len(borders)
        #Список прямоугольников 
        line_list = []
        #Список количества прямоугольников в стороне рамки
        line_counters = []
    
        #Получаем единый массив для всех прямоугольников рамки
        for i in range(borders_length):
            line_counters.append(borders[i].size[0] // gcd_result)
            for j in range(line_counters[i]):
                line_list.append(borders[i].crop(
                        (
                            j * gcd_result, 0, 
                            (j + 1) * gcd_result, border_width
                        )
                    )
                )
    
        #Выполняем циклический сдвиг
        line_list = circular_shift(line_list, step_count)
    
        #Собираем рамку обратно
        for i in range(borders_length):
            for j in range(line_counters[i]):
                borders[i].paste(
                    line_list.pop(0),
                    (j * gcd_result, 0, (j + 1) * gcd_result, border_width)
                )
            for j in range(i):
                borders[i] = borders[i].transpose(Image.ROTATE_90)
            img.paste(borders[i], border_boxes[i])
    
        return img
    
    #----------------------------------------------------------------------------#
    
    img = Image.open(path)
    
    width = img.size[0]
    height = img.size[1]
    
    #Подсчитываем значение ширины рамки
    border_width = int(border_scale * height)
    
    #Коробка для центральной части изображения
    center_box = (
        border_width, border_width,
        width - border_width,  height - border_width
    )
    
    tmp_img = get_new_image(img, width, height, border_width, step_count * (ROW_COUNT - 1))
    #tmp_img.show()
    image_center = tmp_img.crop(center_box)
    tmp_img.paste(image_center.filter(ImageFilter.BLUR), center_box)
    tmp_img.show()
    
def Rand_Rectangle_2(path):
    img = Image.open('photo.png')
    draw = ImageDraw.Draw(img)
    
    
    
    pix = img.load()
    width, height = img.size
    
    #________BOXES______________
    
    draw.line((int((1./6.)*width),0,int((1./6.)*width),height),fill='black',width=1)
    draw.line((int((2./6.)*width),0,int((2./6.)*width),height),fill='black',width=1)
    draw.line((int((3./6.)*width),0,int((3./6.)*width),height),fill='black',width=1)
    draw.line((int((4./6.)*width),0,int((4./6.)*width),height),fill='black',width=1)
    draw.line((int((5./6.)*width),0,int((5./6.)*width),height),fill='black',width=1)
    
    draw.line((0,int((1./6.)*height),width,int((1./6.)*height)),fill='black',width=1)
    draw.line((0,int((2./6.)*height),width,int((2./6.)*height)),fill='black',width=1)
    draw.line((0,int((3./6.)*height),width,int((3./6.)*height)),fill='black',width=1)
    draw.line((0,int((4./6.)*height),width,int((4./6.)*height)),fill='black',width=1)
    draw.line((0,int((5./6.)*height),width,int((5./6.)*height)),fill='black',width=1)
    
    
    
    #NEGATIVE
    for i in range(int((1./6.)*width)):
            for j in range(height):
                a = pix[i, j][0]
                b = pix[i, j][1]
                c = pix[i, j][2]
                draw.point((i, j), (255 - a, 255 - b, 255 - c))
    for i in range(int((2./6.)*width),int((3./6.)*width)):
            for j in range(height):
                a = pix[i, j][0]
                b = pix[i, j][1]
                c = pix[i, j][2]
                draw.point((i, j), (255 - a, 255 - b, 255 - c))
    for i in range(int((4./6.)*width),int((5./6.)*width)):
            for j in range(height):
                a = pix[i, j][0]
                b = pix[i, j][1]
                c = pix[i, j][2]
                draw.point((i, j), (255 - a, 255 - b, 255 - c))
                
                

    depth = 30
    for i in range(int((1./6.)*width),int((2./6.)*width)):
        for j in range(height):
            a = pix[i, j][0]
            b = pix[i, j][1]
            c = pix[i, j][2]
            S = (a + b + c) // 3
            a = S + depth * 2
            b = S + depth
            c = S
            if (a > 255):
                a = 255
            if (b > 255):
                b = 255
            if (c > 255):
                c = 255
            draw.point((i, j), (a, b, c))
    for i in range(int((3./6.)*width),int((4./6.)*width)):
        for j in range(height):
            a = pix[i, j][0]
            b = pix[i, j][1]
            c = pix[i, j][2]
            S = (a + b + c) // 3
            a = S + depth * 2
            b = S + depth
            c = S
            if (a > 255):
                a = 255
            if (b > 255):
                b = 255
            if (c > 255):
                c = 255
            draw.point((i, j), (a, b, c))
    for i in range(int((5./6.)*width),int((6./6.)*width)):
        for j in range(height):
            a = pix[i, j][0]
            b = pix[i, j][1]
            c = pix[i, j][2]
            S = (a + b + c) // 3
            a = S + depth * 2
            b = S + depth
            c = S
            if (a > 255):
                a = 255
            if (b > 255):
                b = 255
            if (c > 255):
                c = 255
            draw.point((i, j), (a, b, c))
    
           
            
    #___________strokes________
    
    line1=(0, 0, width, m.floor((1./6.)*height))
    region1 = img.crop(line1)
        
        #line2=(0,int((1./6.)*height),width,int((2./6.)*height))
        #region2 = Image.crop(line2)
    
    line3=(0,m.ceil((2./6.)*height),width,m.floor((3./6.)*height))
    region3 = img.crop(line3)
        
        #line4=(0,int((3./6.)*height),width,int((4./6.)*height))
        #region4 = Image.crop(line4)
    
    line5=(0,m.ceil((4./6.)*height),width,m.floor((5./6.)*height))
    region5 = img.crop(line5)
    
    
    img.paste(region1,line3)
    img.paste(region3,line5)
    img.paste(region5,line1)
            
            
    
    
    
    
    
    
    
    
    img.show()

def Triangle_Fill(path):
    im1=Image.open(path) # сама картинка 
    x, y = im1.size  # ширина (x) и высота (y) изображения
    img2 = Image.new('RGB', (x, y))
    draw = ImageDraw.Draw(img2)
    
    x1=(0,randint(0,y))
    x3=(randint(0,x),0)
    x5=(x,randint(0,y))
    x7=(randint(0,x),y)
    points1=((x1),(0,0),(x3))
    points2=((x3),(x,0),(x5))
    points3=((x5),(x,y),(x7))
    points4=((x7),(0,y),(x1))
    
    img2.paste(im1, (0,0))
    draw.polygon((points1), fill=(randint(0,255),randint(0,255),randint(0,255)), outline=(0, 0, 0))
    draw.polygon((points2), fill=(randint(0,255),randint(0,255),randint(0,255)), outline=(0, 0, 0))
    draw.polygon((points3), fill=(randint(0,255),randint(0,255),randint(0,255)), outline=(0, 0, 0))
    draw.polygon((points4), fill=(randint(0,255),randint(0,255),randint(0,255)), outline=(0, 0, 0))
    img2.show()

def Rand_Lines(path):
    im = Image.open(path)

    w = im.size[0]
    h = im.size[1]
    
    im1 = im.crop((0, 0, w/10, h))
    im3 = im.crop((w/5, 0, 3*w/10, h))
    im5 = im.crop((2*w/5, 0, w/2, h))
    im7 = im.crop((3*w/5, 0, 7*w/10, h))
    im9 = im.crop((4*w/5, 0, 9*w/10, h))
    
    index = randint(0, 1)
    print (index)
    
    
    if (index == 0):           # нечетные сдвигаются вправо
        im.paste(im1,(w // 10 * 8,0))
        im.paste(im3,(w // 10 * 7,0))
        im.paste(im5,(w // 10 * 6,0))
        im.paste(im7,(w // 10 * 5,0))
        im.paste(im9,(w // 10 * 9,0 ))
    else:                      # нечетные сдвигаются влево
        im.paste(im5,(w // 10,0))
        im.paste(im7,(w // 10 * 2,0))
        im.paste(im9,(w // 10 * 3,0))
        im.paste(im1,(w // 10 * 4,0))
        im.paste(im3,(0,0))
        
    im.show()
    
def Rand_30_rectangles(path):
    im = Image.open(path)

    w = im.size[0]
    h = im.size[1]
    #print (w, h)
    
    m0 = im.crop((0, 0, w/6, h/5)) 
    m1 = im.crop((w/6, 0, w/3, h/5))
    m2 = im.crop((w/3, 0, w/2, h/5))
    m3 = im.crop((w/2, 0, 2*w/3, h/5))
    m4 = im.crop((2*w/3, 0, 5*w/6, h/5))
    m5 = im.crop((5*w/6, 0, w, h/5)) 
    
    
    m6 = im.crop((0, h/5, w/6, 2*h/5))
    
    m7 = im.crop((w/6, h/5, w/3, 2*h/5))
    
    m8 = im.crop((w/3, h/5, w/2, 2*h/5)) 
    m9 = im.crop((w/2, h/5, 2*w/3, 2*h/5))
    
    
    m10 = im.crop((2*w/3, h/5, 5*w/6, 2*h/5))
    
    m11 = im.crop((5*w/6, h/5, w, 2*h/5))
    
    m12 = im.crop((0, 2*h/5, w/6, 3*h/5))
    m13 = im.crop((w/6, 2*h/5, w/3, 3*h/5))
    
    m14 = im.crop((w/3, 2*h/5, w/2, 3*h/5)) 
    
    m15 = im.crop((w/2, 2*h/5, 2*w/3, 3*h/5)) 
    
    m16 = im.crop((2*w/3, 2*h/5, 5*w/6, 3*h/5))
    
    m17 = im.crop((5*w/6, 2*h/5, w, 3*h/5)) 
    
    m18 = im.crop((0, 3*h/5, w/6, 4*h/5))
    
    m19 = im.crop((w/6, 3*h/5, w/3, 4*h/5))
    
    m20 = im.crop((w/3, 3*h/5, w/2, 4*h/5))
    
    m21 = im.crop((w/2, 3*h/5, 2*w/3, 4*h/5))
    
    m22 = im.crop((2*w/3, 3*h/5, 5*w/6, 4*h/5))
    m23 = im.crop((5*w/6, 3*h/5, w, 4*h/5))
    
    
    m24 = im.crop((0, 4*h/5, w/6, h))
    
    m25 = im.crop((w/6, 4*h/5, w/3, h)) 
    
    m26 = im.crop((w/3, 4*h/5, w/2, h))
    
    m27 = im.crop((w/2, 4*h/5, 2*w/3, h)) 
    m28 = im.crop((2*w/3, 4*h/5, 5*w/6, h))
    
    m29 = im.crop((5*w/6, 4*h/5, w, h))
    
    
    im_list = [m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24, m25, m26, m27, m28, m29]
    number_list = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
    n_list = []
    num_list = np.zeros(12)
    i_list = []
    n_list.append(random.sample(number_list, 12))
    for i in range(12):
        a = n_list[0][i]
        i_list.append(im_list[a])
        num_list[i] = n_list[0][i]
    print(num_list)
    random.shuffle(num_list)
    print(num_list)
    
    img1 = Image.new( mode = "RGB", size = (w, h), color = (210, 105, 30) ) 
    
    
    img = img1
    
    def F(image, n):
        if n == 0:
            img.paste(image,(0,0))
        elif n == 1:
            img.paste(image,(int(w/6),0))
        elif n == 2:
            img.paste(image,(int(w/3),0))
        elif n == 3:
            img.paste(image,(int(w/2),0))
        elif n == 4:
            img.paste(image,(int(2*w/3),0))
        elif n == 5:
            img.paste(image,(int(5*w/6),0))
        elif n == 6:
            img.paste(image,(0,int(h/5)))
        elif n == 7:
            img.paste(image,(int(w/6),int(h/5)))
        elif n == 8:
            img.paste(image,(int(w/3),int(h/5)))
        elif n == 9:
            img.paste(image,(int(w/2),int(h/5)))
        elif n == 10:
            img.paste(image,(int(2*w/3),int(h/5)))
        elif n == 11:
            img.paste(image,(int(5*w/6),int(h/5)))
        elif n == 12:
            img.paste(image,(0,int(2*h/5)))
        elif n == 13:
            img.paste(image,(int(w/6),int(2*h/5)))
        elif n == 14:
            img.paste(image,(int(w/3),int(2*h/5)))
        elif n == 15:
            img.paste(image,(int(w/2),int(2*h/5)))
        elif n == 16:
            img.paste(image,(int(2*w/3),int(2*h/5)))
        elif n == 17:
            img.paste(image,(int(5*w/6),int(2*h/5)))
        elif n == 18:
            img.paste(image,(0,int(3*h/5)))
        elif n == 19:
            img.paste(image,(int(w/6),int(3*h/5)))
        elif n == 20:
            img.paste(image,(int(w/3),int(3*h/5)))
        elif n == 21:
            img.paste(image,(int(w/2),int(3*h/5)))
        elif n == 22:
            img.paste(image,(int(2*w/3),int(3*h/5)))
        elif n == 23:
            img.paste(image,(int(5*w/6),int(3*h/5)))
        elif n == 24:
            img.paste(image,(0,int(4*h/5)))
        elif n == 25:
            img.paste(image,(int(w/6),int(4*h/5)))
        elif n == 26:
            img.paste(image,(int(w/3),int(4*h/5)))
        elif n == 27:
            img.paste(image,(int(w/2),int(4*h/5)))
        elif n == 28:
            img.paste(image,(int(2*w/3),int(4*h/5)))
        elif n == 29:
            img.paste(image,(int(5*w/6),int(4*h/5)))
        
    for i in range(12):
        image = i_list[i]
        n = num_list[i]
        F(image, n)
    img.show()
    
def TR_2(path):
    pic1 = Image.open(path)
    objVP=pic1.load()
    w=853 
    h=853 # длина и ширина фотки
    n=0
    m=0
    k=50 #коэффициент сепии
    for j in range(0,283,1):
        for i in range(0,283,1):
            n=284-i-j
            m=852-284+i+j
            if (m<0 or n<0): break 
            r=objVP[284-i-j,852-j][0]
            g=objVP[284-i-j,852-j][1]
            b=objVP[284-i-j,852-j][2]
            S=(r+g+b)//3
            r1=S+2*k
            g1=S+k
            b1=S
            if (r1 > 255): r1 = 255
            if (g1 > 255): g1 = 255
            if (b1 > 255): b1 = 255 
            objVP[852-284+i+j,j]=(r1,g1,b1)
    pic1Crop=pic1.crop((0,0,w,400))  
    
    m=0
    n=0
    pic2 = Image.open(path)     
    objNL=pic2.load()
    for j in range(0,283,1):
        for i in range(0,283,1):
            n=284-i-j
            m=852-284+i+j
            if (m<0 or n<0): break 
            r=objNL[852-284+i+j,j][0]
            g=objNL[852-284+i+j,j][1]
            b=objNL[852-284+i+j,j][2]
            S=(r+g+b)//3
            r1=S+2*k
            g1=S+k
            b1=S
            if (r1 > 255): r1 = 255
            if (g1 > 255): g1 = 255
            if (b1 > 255): b1 = 255
            objNL[284-i-j,852-j]=(r1,g1,b1)
    pic2Crop=pic2.crop((0,401,w,853))
    
    pic3=Image.new('RGB',(w,h))
    pic3.paste(pic1Crop,(0,0))
    pic3.paste(pic2Crop,(0,400))
    # ДЕЛИМ ИЗОБРАЖЕНИЕ НА 3 ЧАСТИ ДЛЯ УДОБСТВА И ПРИМЕНЯЕМ ЧБ ОТДЕЛЬНО
    obj1=pic3.load()
    for j in range(0,283,1):
        for i in range(0,852-283+j,1):
            r=obj1[i,j][0]
            g=obj1[i,j][1]
            b=obj1[i,j][2]
            S=(r+g+b)//3
            obj1[i,j]=(S,S,S)
              
    for j in range(283,852-283,1):
        for i in range(0,852,1):
            r=obj1[i,j][0]
            g=obj1[i,j][1]
            b=obj1[i,j][2]
            S=(r+g+b)//3
            obj1[i,j]=(S,S,S)           
    
    for j in range(852-284,852,1):
        for i in range(-852+284+j,852,1):
            r=obj1[i,j][0]
            g=obj1[i,j][1]
            b=obj1[i,j][2]
            S=(r+g+b)//3
            obj1[i,j]=(S,S,S)            
    
    
        
    pic3.show()  

def TR_3(path):
    

    img = Image.open(path) #объявляем картинку
    img.save("save.png", "PNG")
    img = Image.open('save.png')
    img = img.convert("RGBA")
    newData_2 = []
    img.size
    pixdata = img.load()
    
    width_prob=0
    height_prob=0
    
    width_prob=img.size[0]
    height_prob=img.size[1]
    
    shirina=width_prob
    visota=height_prob
    
    with Image.open("save.png") as img:
        if (width_prob<height_prob):
            razn=height_prob-width_prob
            (left, upper, right, lower) = (0, 0, width_prob, height_prob)
            img_crop_3 = img.crop((left, upper, right+razn, lower))
            draw = ImageDraw.Draw(img_crop_3) 
            draw.polygon([(width_prob,0), (width_prob+razn,0), (width_prob+razn,height_prob), (width_prob, height_prob)], fill = '#ffffff')
            img_crop_3 = img_crop_3.rotate(90)
            (left, upper, right, lower) = (0, razn, height_prob, razn+width_prob)
            img_crop_3 = img_crop_3.crop((left, upper, right, lower))
            img_crop_3.save('save.png')
    
    img = Image.open('save.png')
    img = img.convert("RGBA")
    newData_2 = []
    img.size
    pixdata = img.load()
    
    width=img.size[0]
    height=img.size[1]
    
    
    b=(height*height)/(2*width)
    x=height-(width/2)
    
    
    with Image.open("save.png") as img: #рисуем белые треугольники
        #2
        #обрезка
        (left, upper, right, lower) = (width/2, 0, width, height)
        img_crop_2 = img.crop((left, upper, right+x, lower))
        draw = ImageDraw.Draw(img_crop_2) 
        
        #треугольники
        draw.polygon([(0,0), (width/2, 0), (0,height/2)], fill = '#ffffff')
        draw.polygon([(0,height), (width/2, height), (0,height/2)], fill = '#ffffff')
        draw.polygon([(width/2,0), ((width/2)+x,0), ((width/2)+x,height), (width/2, height)], fill = '#ffffff')
        img_crop_2.save('tr_2.png')
        
        #1
        (left, upper, right, lower) = (0, 0, width/2, height)
        img_crop_1 = img.crop((left, upper, right+x, lower))
        draw = ImageDraw.Draw(img_crop_1)
        draw.polygon([(0,0), (width/2, 0), (width/2,height/2)], fill = '#ffffff')
        draw.polygon([(0,height), (width/2, height), (width/2,height/2)], fill = '#ffffff')
        draw.polygon([(width/2,0), ((width/2)+x,0), ((width/2)+x,height), (width/2, height)], fill = '#ffffff')
        img_crop_1.save('tr_1.png')
    
    tr_2 = Image.open('tr_2.png')
    rotated_img_2 = tr_2.rotate(-90)
    (left, upper, right, lower) = (0, 0, height, (height)/2)
    rotated_img_2 = rotated_img_2.crop((left, upper, right, lower))
    
    draw=ImageDraw.Draw(rotated_img_2)
    pixelm=rotated_img_2.load()
    
    for i in range(int(height)): #картинка чб
        for j in range(int(height/2)):
            a=pixelm[i,j][0]
            b=pixelm[i,j][1]
            c=pixelm[i,j][2]
            s1=(a+b+c)//3
            draw.point((i,j),(s1,s1,s1))
    
    rotated_img_2.save("tr_2.png", "PNG")
    
    tr_1 = Image.open('tr_1.png')
    rotated_img_1 = tr_1.rotate(-90)
    (left, upper, right, lower) = (0, height/2-x, height, (height)-x)
    rotated_img_1 = rotated_img_1.crop((left, upper, right, lower))
    #rotated_img_1.save("tr_1.png", "PNG")
    
    draw=ImageDraw.Draw(rotated_img_1)
    pixelm=rotated_img_1.load()
    
    for i in range(int(height)): #картинка негатив
        for j in range(int(height/2)):
            a=pixelm[i,j][0]
            b=pixelm[i,j][1]
            c=pixelm[i,j][2]
            draw.point((i,j),(255-a, 255-b, 255-c))
            
    rotated_img_1.save("tr_1.png", "PNG")
    
    #прозрачный
    tr_1 = Image.open('tr_1.png')
    tr_1 = tr_1.convert("RGBA")
    datas_1 = tr_1.getdata()
    newData_1 = []
    for item in datas_1:
        if item[0] == 0 and item[1] == 0 and item[2] == 0:
            newData_1.append((0, 0, 0, 0))
        else:
            newData_1.append(item)
    tr_1.putdata(newData_1)
    
    tr_1.save("tr_1.png", "PNG")
    
    
    tr_2 = Image.open('tr_2.png')
    tr_2 = tr_2.convert("RGBA")
    datas_2 = tr_2.getdata()
    newData_2 = []
    for item in datas_2:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData_2.append((0, 0, 0, 0))
        else:
            newData_2.append(item)
    tr_2.putdata(newData_2)
    tr_2.save("tr_2.png", "PNG")
    
    
    img = Image.open('save.png')
    img = img.convert("RGBA")
    
    tr_1 = Image.open('tr_1.png')
    tr_1 = tr_1.convert("RGBA")
    
    img.paste(tr_1, (int((width-height)/2), 0),  tr_1)  
    
    tr_2 = Image.open('tr_2.png')
    tr_2 = tr_2.convert("RGBA")
    
    img.paste(tr_2, (int((width-height)/2), int(height/2)),  tr_2)  
        
    
    img.save("itog.png", "PNG")
    
    img = Image.open('itog.png')
    if (shirina<visota):
        
        img = img.convert("RGBA")
        newData_2 = []
        img.size
        pixdata = img.load()
    
        width_prob=img.size[0]
        height_prob=img.size[1]
    
    
        with Image.open("itog.png") as img:
            razn=abs(height_prob-width_prob)
            (left, upper, right, lower) = (0, 0, width, height)
            img_crop_3 = img.crop((left, upper, right, lower+razn))
            #draw = ImageDraw.Draw(img_crop_3) 
            #draw.polygon([(width_prob,0), (width_prob+razn,0), (width_prob+razn,height_prob), (width_prob, height_prob)], fill = '#ffffff')
            img_crop_3 = img_crop_3.rotate(-90)
            (left, upper, right, lower) = (razn, 0, razn+height, width)
            img_crop_3 = img_crop_3.crop((left, upper, right, lower))
            img_crop_3.save('itog.png')
            img_crop_3.show()
    
    img.show()
    os.remove('itog.png')
    os.remove('tr_2.png')
    os.remove('tr_1.png')
    os.remove('save.png')

def Donut(path):
    img = Image.open(path) #объявляем картинку
    img.save("shenki.png", "PNG")
    img = Image.open('shenki.png')
    img = img.convert("RGBA")
    newData_2 = []
    img.size
    pixdata = img.load()
    width=img.size[0]
    height=img.size[1]
    
    
    im_maska = Image.new('RGB', (width, height), (32, 178, 170))
    
    
    #20B2AA
    
    if(width > height):
        r_out = height/2
    else:
        r_out = width/2
    
    r_in = int(random.uniform((1/3)*r_out, (1/4)*r_out))
    
    
    draw = ImageDraw.Draw(im_maska)
    draw.ellipse((int((width/2)-r_out), 0, int((width/2)+r_out), height), fill='#ffffff')
    draw.ellipse((int((width/2)-r_in), int(height/2)-r_in, int((width/2)+r_in), int(height/2)+r_in), fill='#20B2AA')
    im_maska.save('maska.png')
    
    prozr = Image.open('maska.png')
    prozr = prozr.convert("RGBA")
    datas = prozr.getdata()
    
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    
    prozr.putdata(newData)
    prozr.save("maska.png", "PNG")
    
    img = Image.open('shenki.png')
    img = img.convert("RGBA")
    
    prozr = Image.open('maska.png')
    prozr = prozr.convert("RGBA")
    
    img.paste(prozr, (0, 0),  prozr)
    
    img.show()
    
    os.remove('maska.png')
    os.remove('shenki.png')
    
def Four_Rectangles(path):
    im = Image.open(path)
    w, h = im.size
    fon = Image.new( mode = "RGB", size = (w, h) )
    
    im1=im.crop((0,0,w // 3,h // 2)) #обрезка
    im2=im.crop((w // 3,h // 2,w,h))
    im3=im.crop((w // 3, 0, w, h // 2))
    im4=im.crop((0, h // 2,w // 3, h))
    
    
    draw=ImageDraw.Draw(im2)
    width=im2.size[0]
    height=im2.size[1]
    pix=im2.load()
    
    for i in range(width):
        for j in range(height):
            a=pix[i,j][0]
            b=pix[i,j][1]
            c=pix[i,j][2]
            draw.point((i,j),(255-a,255-b,255-c))
            
    
    draw=ImageDraw.Draw(im3)
    width=im3.size[0]
    height=im3.size[1]
    pix=im3.load()
    depth=30
    
    for i in range(width):
        for j in range(height):
            a=pix[i,j][0]
            b=pix[i,j][1]
            c=pix[i,j][2]
            s=(a+b+c)//3
            a=s+depth*2
            b=s+depth
            c=s
            if(a>255):
                a=255
            if(b>255):
                b=255
            if(c>255):
                c=255
            draw.point((i,j),(a,b,c))
            
    fon.paste(im1,(w - w // 3, h // 2))
    
    pc1 = fon
    pc1.paste(im2,(0,0))
    
    
    pc1.paste(im4,(w - w // 3,0))
    
    
    pc1.paste(im3,(0,h // 2))
    pc1.show()

def Heart(path):

    def cardioid_coordinates(r: int, t: float, w, h) -> (float, float):
        x = 2 * r * m.sin(t) + r * m.sin(t * 2) + w / 2 
        y = 2 * r * m.cos(t) + r * m.cos(t * 2) + h / 2.5
        return x, y


    def cardioid(r, step, w, h):
        coordinates = [cardioid_coordinates(r, i/step * m.pi, w, h) for i in range(0, 2 * step, 1)]
        img = Image.new('RGBA', (w, h), 'pink')
        cardioid_draw = ImageDraw.Draw(img)
        cardioid_draw.polygon(coordinates, fill='white')
        return img
    pc4 = Image.open('photo.png')
    w,h = pc4.size
    img = cardioid(w // 7,h // 7, w, h)
    
    frm = img
    img = frm.convert("RGBA")
    datas = img.getdata()
    #делает все белые пиксели прозрачными
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    
    frm.putdata(newData)
    
    
    pc4.paste(frm,mask=frm)
    
    pc4.show()
