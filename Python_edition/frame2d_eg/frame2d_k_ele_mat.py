import numpy as np
from math import *

def frame2d_k_ele_mat(ele_info):
    #生成局部坐标系下的单元刚度矩阵，ele_info形式为[e;a;iz;l]。
    e = ele_info[0]
    a = ele_info[1]
    iz = ele_info[2]
    l = ele_info[3]
    x1 = e * a / l
    x2 = 12.0 * e * iz / pow(l,3)
    x3 = 6.0 * e * iz / pow(l,2)
    x4 = 4.0 * e * iz / l
    x5 = 2.0 * e * iz / l
    aa = np.array([[x1, 0, 0, -x1, 0, 0],
    [0, x2, x3, 0, -x2, x3],
    [0, x3, x4, 0, -x3, x5],
    [-x1, 0, 0, x1, 0, 0],
    [0, -x2, -x3, 0, x2, -x3],
    [0, x3, x5, 0, -x3, x4]])
    y = aa

    return y