import numpy as np
from math import *
import sympy as sp

def frame2d_work_equi_distributed(l,q):

    #分布力载荷在节点产生的等效力, 此处特殊考虑均匀横向分布力，大小为q。
    x=sp.symbols('x')
    n = np.array([1 - 3 * x**2 / l**2 + 2 * x** 3 / l** 3,
    x - 2 * x** 2 / l + x **3 / l** 2,
    3 * x **2 / l** 2 - 2 * x **3 / l** 3,
    x **3 / l ** 2 - x ** 2 / l])
    f = np.zeros((4, 1))
    for i in range(4):
        f[i] = sp.integrate(n[i] * q, (x,0, l))

    return f
