# Этот код имитирует осциллятор затухания:
# Осциллятор с затухающими гармониками в двухъямном потенциале.

# F = -gamma*dx/dt + 2*a*x - 4*b*x^3 + F_0*cos(omega*t)

# Нелинейное дифференциальное уравнение второго порядка, решенное численно с помощью разложения Тейлора.

# Для текущего набора параметров движение хаотическое, т.е.
# движение сильно зависит от начальных условий. 
# Дополнительно фиксированного периода движения не наблюдается. Сюжет Пуанкаре - фрактал.

from numpy import * 
from matplotlib.pyplot import * 

# -------------Коэффициенты---------------------

# Первоначальные условия
x,v = 0.5, 0

# Продолжительность моделирования
T = 1000

# Для хаоса а = 0.5 и b = 0.0625
# Для простого случая а = 0 и b = 1

a = 1	
b = 1 	

F_0 = 2.5

omega = 2.0
gamma = 0.1

# ------------------------------------------

h = 1e-1 # Шаг времени

period = 2*pi/(1.0*omega)
t = arange(0,T,h)


def x_2(x,v):
    return -gamma*v + 2.0*a*x - 4.0*b*x*x*x

def x_3(x2,x,v):
    return -gamma*x2 + 2.0*a*v -12.0*b*x*x*v

def x_4(x3,x2,x,v):
    return -gamma*x3 + 2.0*a*x2 -12.0*b*x*x*x2 - 24.0*b*v*v*x

def x_5(x4,x3,x2,x,v):
    return -gamma*x4 + 2*a*x3 -12.0*b*(x*x*x3 + 2.0*x2*x*v) -24.0*b*(v*v*v+2*x*v*x2)

# Тригонометрические члены в производных. Оценить перед циклом
x2F = F_0*cos(omega*t)
x3F = -F_0*omega*sin(omega*t)
x4F = -F_0*omega*omega*cos(omega*t)
x5F = F_0*omega*omega*omega*sin(omega*t)

# Коэффициенты перед разложением в ряд Тейлора
coef1 = 0.5*h**2.0
coef2 = 1.0/6.0*h**3.0
coef3 = 1.0/24.0*h**4.0
coef4 = 1.0/120.0*h**5.0


position, velocity = zeros(len(t)) , zeros(len(t))
position[0] = x

for i in range(1,len(t)):

    d2 = x_2(x,v) + x2F[i]
    d3 = x_3(d2,x,v) + x3F[i]
    d4 = x_4(d3,d2,x,v) + x4F[i]
    d5 = x_5(d4,d3,d2,x,v) + x5F[i]
	
    # Разложение в ряд Тейлора для x, v. Заказать h ^ 5
    x += v*h + coef1*d2 + coef2*d3 + coef3*d4 + coef4*d5
    v += d2*h + coef1*d3 + coef2*d4 + coef3*d5
	
    position[i], velocity[i] = x , v


# Получите точки фазового пространства в целых кратных периодах для графика Пуанкаре
strange_attractor = zeros([int(T/period),2])
k = 1

for i in range(len(t)):
    if abs(t[i]-k*period)<h:
        strange_attractor[k-1,0] = position[i]
        strange_attractor[k-1,1] = velocity[i]
        k+=1


subplot(3,1,1)
plot(t[-3000:],position[-3000:],'g-')
title('Траектория осциллятора')
ylabel('Позиция')

subplot(3,1,2)
plot(position[-3000:],velocity[-3000:],'r-')
title('Фазовое пространство')
xlim([-4.5,4.5])
ylabel('Момент')

subplot(3,1,3)
scatter(strange_attractor[:,0],strange_attractor[:,1])
xlabel('Позиция')
ylabel('Момент')
title(r'График Пуанкаре (Фазовое пространство во времени = $\frac{2\pi N}{\omega}$, N = 1,2,3...)')

show()
