import numpy as np
from math import *
import matplotlib.pyplot as plt
import sympy as sp

def frame2d_beam_deflection_plot(l, u, eles, ele,axes_def,axes_slo):

    #绘出各单元的挠度曲线，u为四行一列。
    x=sp.symbols('x')
    deflection=1000*np.dot(np.array([1-3*x**2/l**2+2*x**3/l**3,x-2*x**2/l+x**3/l**2,3*x**2/l**2-2*x**3/l**3,x**3/l**2-x**2/l]),u)
    slope=sp.diff(deflection,x)/1000*180/pi
    length=np.linspace(0,l,l/0.01+1)
    b=[deflection.subs(x,length[i]) for i in range(length.shape[0])]
    c=[slope.subs(x,length[i]) for i in range(length.shape[0])]

    axes_def.plot(length,b)
    axes_def.set_ylabel('Deflection/mm')
    axes_def.grid()
    axes_slo.plot(length,c)
    axes_slo.set_ylabel('slope/degree')
    axes_slo.grid()


