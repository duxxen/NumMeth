import math
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from typing import Callable

''' 0) Исходные данные '''

C = -30
N = 9
K = 83

a = np.pi * (N + 10)
b = a + K / 50 + 2

h = 0.0001
n = int((b - a) / h)
h = (b - a) / (n)

''' 1) Параметры: p(x), q(x) '''

def p(x: float):
    return np.sin((N * x) / np.pi) + x * np.cos(2 * N) + C

def dp(x: float):
    return np.cos(N * x / np.pi) * N / np.pi + np.cos(2 * N)

def q(x: float):
    return np.exp(np.cos(N * x / np.pi)) + np.pow(np.cos(N * x / 10), 2) + 1

''' 2) Нетривиальный физически осмысленный пример: u(x), f(x) '''

def u(x: float):
    return np.pow(np.sin(np.pi * (x - b) / (b - a)), 2)

def du(x: float):
    return np.pi * np.sin(2 * np.pi * (x - b) / (b - a)) / (b - a)

def ddu(x: float):
    return 2 * np.pow(np.pi, 2) * np.cos(2 * np.pi * (x - b) / (b - a)) / np.pow(b - a, 2)

def f(x: float):
    return -dp(x) * du(x) - p(x) * ddu(x) + q(x) * u(x)

''' 3)  Метод МКЭ '''
''' 3.1) Методы рассчета интеграла '''

def quad(f: Callable[[float], float], a: float, b: float, h: float):
    return integrate.quad(f, a, b)[0]

def simp(f: Callable[[float], float], a: float, b: float, h: float):
    n = int((b - a) / h)
    h = (b - a) / n
    x = np.linspace(a, b, n)

    s1 = 0
    s2 = 0
    for i in range(2, n-1, 2):
        s1 += 4 * f(x[i])
    for i in range(1, n-1, 2):
        s2 += 2 * f(x[i])

    return h/3 * (f(x[0]) + s1 + s2 + f(x[-1]))

''' 3.2) Конечные функции phi(x) '''

def xk(k: int) -> float:
    return max(min(a + k * h, b), a)

def phi(x: float, k: int) -> float:
    hn = h * np.sqrt(h)
    if xk(k-1) <= x and x <= xk(k):     return (x - xk(k-1)) / hn
    elif xk(k) < x and x <= xk(k+1):    return (xk(k+1) - x) / hn
    else:                               return 0

''' 3.3) Генерация СЛАУ '''

def createMatrix(ifunc = quad):
    mtr = np.zeros((n-2, n-2))
    md = np.zeros(n-2)        # главная i
    bd = np.zeros(n-3)      # под   i-1
    ad = np.zeros(n-3)      # над   i+1

    h3 = 1 / np.pow(h, 3)

    for k in range(0, n-2):
        xkp = xk(k)
        xki = xk(k+1)
        xkn = xk(k+2)
        md[k] = ifunc(lambda x: -h3 * p(x) + q(x) * np.pow(phi(x, k+1), 2), xkp, xkn, h/10)

        if k == n-3: break
        
        bd[k] = ifunc(lambda x: -h3 * p(x) + q(x) * phi(x, k+1) * phi(x, xki), xkp, xki, h/10)
        ad[k] = ifunc(lambda x: -h3 * p(x) + q(x) * phi(x, k+2) * phi(x, xkn), xki, xkn, h/10)
    
    np.fill_diagonal(mtr, md)
    np.fill_diagonal(mtr[:, 1:], ad)
    np.fill_diagonal(mtr[1:, :], bd)

    return mtr

def createVect(ifunc = quad):
    vct = np.zeros(n-2)
    for k in range(n-2):
        vct[k] = ifunc(lambda x: f(x) * phi(x, k+1), xk(k), xk(k+2), h/10)
    
    return 1/np.sqrt(h) * vct

''' 3.4) Метод решения СЛАУ из ЛР1 - Метод Зейделя '''

def SeidelMethod(mtr: np.ndarray[float], vct: np.ndarray[float], x0: np.ndarray[float], eps = 1e-10):
    n = len(mtr)
    xp = x0.copy()
    converge = False

    while not converge:
        x = xp.copy()
        for i in range(n):
            s1 = 0
            s2 = 0
            if (i > 0):         s1 = mtr[i][i - 1] * x[i - 1]
            elif (i < n - 1):   s2 = mtr[i][i + 1] * xp[i + 1]

            x[i] = (vct[i] - s1 - s2) / mtr[i][i]
        xp = x
        converge = np.linalg.norm(x - xp) <= eps
    return x

''' 3.5) Функция выполнения ЛР '''

def lab2(type: str, ifunc = Callable[[float], float]):

    global n, h
    fig, (axtst, axrh) = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle(f'Решение методом {type}')
    axtst.set_title('Функция u(x)')
    axrh.set_title(f'Погрешность r(h)')

    ''' 4) Проверка работоспособности метода '''

    xv = np.linspace(a, b, n) 
    uv = u(xv)

    matr = createMatrix(ifunc)
    vect = createVect(ifunc)
    uh = [u(a)] + list(np.linalg.solve(matr, vect)) + [u(b)]
    r = np.sqrt(h) * np.linalg.norm(uv - uh)
    print(f'r(h={h}) = {r}')

    axtst.plot(xv, uv, color='red', label='u(x)')
    axtst.plot(xv, uh, color='blue', label=f'u(h={h})')
    axtst.grid()
    axtst.legend()

    ''' 5) Многократное решение системы с изменением h = 10, 20, ..., 1000 '''

    t = 5
    hv = np.zeros(t)
    rv = np.zeros(t)
    for i in range(t):
        h = 1 / (10 * (i + 1))
        n = int((b - a) / h)

        xv = np.linspace(a, b, n)
        uv = u(xv)

        matr = createMatrix(ifunc)
        vect = createVect(ifunc)
        uh = [u(a)] + list(SeidelMethod(matr, vect, vect)) + [u(b)]

        r = np.sqrt(h) * np.linalg.norm(uv - uh)
        hv[i] = h
        rv[i] = r

        print(f'[{type}]: Computing... [{i}/{t}], n={n}')
    
    ''' 6) Построение графика погрешности r(h) '''

    axrh.plot(hv, rv, color='red')
    axrh.set_yscale('log')

    ''' 7) Определение порядка сходимости '''
    print(f'h[0] = {hv[0]}\nh[1] = {hv[1]}')
    print(f'r[0] = {rv[0]}\nr[1] = {rv[1]}')
    print(f'p = {np.sqrt(rv[1]/rv[0])}')

    plt.show()
    


if __name__ == '__main__':

    lab2('quad()', quad)

