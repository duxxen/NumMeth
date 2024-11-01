import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from typing import Callable, overload
from numpy.typing import NDArray
from multipledispatch import dispatch

N = 9
K = 83
C = -30
a = np.pi * (N + 10)
b = a + K / 10 + 2

pc = K * np.exp(10*N / K)
qc = N * np.sin(K ** N) + 2*K

def u(x): return np.sin(np.pi * (x - b) / (b - a)) ** 2
def du(x): return  np.pi * np.sin(2 * np.pi * (-x + b) / (b - a)) / (b - a)
def ddu(x): return 2 * np.pow(np.pi, 2) * np.cos(2 * np.pi * (x - b) / (b - a)) / np.pow(b - a, 2)

def p(x): return np.sin((N * x) / np.pi) + x * np.cos(2 * N) + C
def dp(x): return np.cos((N * x) / np.pi) * N / np.pi
def q(x): return np.exp(np.cos((N * x) / np.pi)) + np.pow(np.cos((N * x) / 10), 2) + 1

def fc(x): return -pc * ddu(x) + qc * u(x)
def f(x): return -dp(x) * du(x) - p(x) * ddu(x) + q(x) * u(x)


def quad(f: Callable[[float], float], a: float, b: float, h: float):
    return integrate.quad(f, a, b)[0]


def simp(f: Callable[[float], float], a: float, b: float, h: float):
    n = int((b - a) / h)
    h = (b - a) / n
    x = np.linspace(a, b, n)
    return h / 3 * f(a) + 4 * f(x[1:-1][::2]).sum() + 2 * f(x[1:-1][1::2]).sum() + f(b)


def createRightVect(f: Callable[[float], float], n: int, ifunc = quad) -> NDArray[np.float64]:
    h = (b - a) / n
    vb = np.zeros(n)

    for k in range(0, n):
        xk = a + (k + 1) * h
        xp = xk - h
        xn = xk + h
        vb[k] = ifunc(lambda x: f(x) * (x - xp), xp, xk, h/10) + ifunc(lambda x: f(x) * (xn - x), xk, xn, h/10)

    return 1/(h**2) * vb


def createMatrixConst(n: int) -> NDArray[np.float64]:
    h = (b - a) / n

    ma1 = np.zeros((n, n))
    np.fill_diagonal(ma1, 2)
    np.fill_diagonal(ma1[1:, :], -1)
    np.fill_diagonal(ma1[:, 1:], -1)
    ma1 *= pc / (h**2)

    ma2 = np.zeros((n, n))
    np.fill_diagonal(ma2, 4)
    np.fill_diagonal(ma2[1:, :], 1)
    np.fill_diagonal(ma2[:, 1:], 1)
    ma2 *= qc / 6

    return ma1 + ma2


def createMatrixFunc(n: int, ifunc = quad) -> NDArray[np.float64]:
    h = (b - a) / n
    ma = np.zeros((n, n))
    mc = 1/(h**2)           # matrix coef
    md = np.zeros(n)        # main diagonal
    ad = np.zeros(n - 1)    # above main
    ud = np.zeros(n - 1)    # under main
    for k in range(0, n):
        xk = a + (k + 1) * h
        xp = xk - h
        xn = xk + h

        md = mc * (ifunc(lambda x: p(x), xp, xn, h/10) + ifunc(lambda x: q(x) * pow(x - xp, 2), xp, xk, h/10) + ifunc(lambda x: q(x) * pow(xn - x, 2), xk, xn, h/10))
        if k == n-1: break
        ad = mc * (-ifunc(lambda x: p(x), xk, xn, h/10) + ifunc(lambda x: q(x) * (x - xk) * (xn - x), xk, xn, h/10))
        ud = mc * (-ifunc(lambda x: p(x), xp, xk, h/10) + ifunc(lambda x: q(x) * (x - xp) * (xk - x), xp, xk, h/10))
    
    np.fill_diagonal(ma, md)
    np.fill_diagonal(ma[:, 1:], ad)
    np.fill_diagonal(ma[1:, :], ud)
    return ma


def SeidelMethod(a: NDArray[np.float64], b: NDArray[np.float64], x0: NDArray[np.float64], eps: float) -> NDArray[np.float64]:
    n = len(a)
    xp = x0.copy()
    converge = False

    while not converge:
        x = xp.copy()
        for i in range(0, n):
            s1 = 0
            s2 = 0
            if (i > 0):         s1 = a[i][i - 1] * x[i - 1]
            elif (i < n - 1):   s2 = a[i][i + 1] * xp[i + 1]

            x[i] = (b[i] - s1 - s2) / a[i][i]
        xp = x
        converge = np.linalg.norm(x - xp) <= eps
    return x

inputFig, inputAxes = plt.subplots(2, 2, figsize=(10, 8))
inputFig.suptitle('Input data')

constTestFig, constTestAxes = plt.subplots(2, figsize=(10, 8))
constTestFig.suptitle('Output data (p,q = const)')

funcTestFig, funcTestAxes = plt.subplots(2, 2, figsize=(10, 8))
funcTestFig.suptitle('Output data (p,q = func)')


def lab2(type: str, h: float, ifunc: Callable[[float], float]):
    n = int((b - a) / h)
    xv = np.linspace(a, b, n)
    uv = u(xv)
    t = 100

    if type is 'const':
        matr = createMatrixConst(n)
        vect = createRightVect(fc, n)
        ux = np.linalg.solve(matr, vect)

        constTestAxes[0].set_title('Function u(x)')
        constTestAxes[0].plot(x, uv, color='red', label='u(x)')
        constTestAxes[0].plot(x, ux, color='blue', label='u(h = {})'.format(h))
        constTestAxes[0].legend()
        constTestAxes[0].grid()

        hv = np.zeros(t)
        rv = np.zeros(t)
        for i in range(0, t):
            h = 1 / (10 * (i + 1))
            n = int((b - a) / h)

            xv = np.linspace(a, b, n)
            uv = u(xv)

            matr = createMatrixConst(n)
            vect = createRightVect(fc, n)
            ux = SeidelMethod(matr, vect, vect, 1e-10)

            r = np.sqrt(h) * np.linalg.norm(uv - ux)
            hv[i] = h
            rv[i] = r

            print('[const]: Computing... [{}/{}], n={}'.format(i+1, t, n))
        
        constTestAxes[1].set_title('Error r(h)')
        constTestAxes[1].plot(hv, rv, color='red', label='r(h)')
        constTestAxes[1].set_xlim(hv.max(), hv.min())
        constTestAxes[1].set_xscale('log')
        constTestAxes[1].minorticks_off()
        constTestAxes[1].grid()
    else:
        method = 'quad'
        col = 0
        if ifunc is simp:
            method = 'Simpson'
            col = 1

        matr = createMatrixFunc(n, ifunc)
        vect = createRightVect(f, n, ifunc)
        ux = np.linalg.solve(matr, vect)

        funcTestAxes[0, col].set_title('Function u(x) ({})'.format(method))
        funcTestAxes[0, col].plot(x, uv, color='red', label='u(x)')
        funcTestAxes[0, col].plot(x, ux, color='blue', label='u(h = {})'.format(h))
        funcTestAxes[0, col].legend()
        funcTestAxes[0, col].grid()

        hv = np.zeros(t)
        rv = np.zeros(t)
        for i in range(0, t):
            h = 1 / (10 * (i + 1))
            n = int((b - a) / h)

            xv = np.linspace(a, b, n)
            uv = u(xv)

            matr = createMatrixFunc(n, ifunc)
            vect = createRightVect(f, n, ifunc)
            ux = SeidelMethod(matr, vect, vect, 1e-10)

            r = np.sqrt(h) * np.linalg.norm(uv - ux)
            hv[i] = h
            rv[i] = r

            print('[{}]: Computing... [{}/{}], n={}'.format(method, i+1, t, n))
        
        funcTestAxes[1, col].set_title('Error r(h) ({})'.format(method))
        funcTestAxes[1, col].plot(hv, rv, color='red', label='r(h)')
        funcTestAxes[1, col].set_xlim(hv.max(), hv.min())
        funcTestAxes[1, col].set_xscale('log')
        funcTestAxes[1, col].minorticks_off()
        funcTestAxes[1, col].grid()




if __name__ == '__main__':
    hbase = 1 / 100
    n = int((b - a) / hbase)

    x = np.linspace(a, b, n)

    inputAxes[0, 0].set_title('Function f(x) (p,q = const)')
    inputAxes[0, 0].plot(x, fc(x), color='red', label='f(x)')
    inputAxes[0, 0].grid()

    inputAxes[0, 1].set_title('Function f(x) (p,q = funcs)')
    inputAxes[0, 1].plot(x, f(x), color='blue', label='f(x)')
    inputAxes[0, 1].grid()

    inputAxes[1, 0].set_title('Function p(x)')
    inputAxes[1, 0].plot(x, p(x), color='orange', label='p(x)')
    inputAxes[1, 0].plot(x, dp(x), label="p'(x)")
    inputAxes[1, 0].legend()
    inputAxes[1, 0].grid()

    inputAxes[1, 1].set_title('Function q(x)')
    inputAxes[1, 1].plot(x, q(x), color='magenta')
    inputAxes[1, 1].grid()

    lab2('const', hbase, None)
    lab2('funcs', hbase, quad)
    lab2('funcs', hbase, simp)

    
    plt.show()

    