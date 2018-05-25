import numpy as np
from math import *


def frame3d_r_mat(ele_info):
    # 生成各单元的转移矩阵。取局部坐标系的原则是y轴在总体坐标系XOY平面内.ele_info[c_Xx;c_Yx;c_Zx]。
    # 一种局部坐标系的取法原则是局部x与全局Z叉乘得局部y坐标轴方向，即本函数里的取法。
    cxx = ele_info[0]
    cyx = ele_info[1]
    czx = ele_info[2]
    if cxx == 0 and cyx == 0:
        if czx > 0:
            r = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
        else:
            r = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    else:
        d = sqrt(cxx ** 2 + cyx ** 2)
        cxy = -cyx / d
        cyy = cxx / d
        czy = 0
        cxz = -czx * cyy
        cyz = -czx * cyx / d
        czz = d
        r = np.array([[cxx, cyx, czx], [cxy, cyy, czy], [cxz, cyz, czz]])

    r0 = np.zeros((3, 3))
    rr1 = np.column_stack((r, r0, r0, r0))
    rr2 = np.column_stack((r0, r, r0, r0))
    rr3 = np.column_stack((r0, r0, r, r0))
    rr4 = np.column_stack((r0, r0, r0, r))
    rr = np.row_stack((rr1, rr2, rr3, rr4))
    y = rr

    return y
