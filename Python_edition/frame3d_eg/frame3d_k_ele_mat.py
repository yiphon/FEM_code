import numpy as np
from math import *

def frame3d_k_ele_mat(ele_info):
    #完成局部坐标系下的单位刚度矩阵，ele_info为主函数中构造的完整的单元信息矩阵。
    e = ele_info[0]
    a = ele_info[1]
    iy = ele_info[2]
    iz = ele_info[3]
    j = ele_info[4]
    l = ele_info[5]
    g = ele_info[9]
    x1 = e * a / l
    x2 = 12 * e * iz / l** 3
    x3 = 6 * e * iz / l ** 2
    x4 = 12 * e * iy / l ** 3
    x5 = 6 * e * iy / l ** 2
    x6 = g * j / l
    x7 = 4 * e * iy / l
    x8 = 2 * e * iy / l
    x9 = 4 * e * iz / l
    x10 = 2 * e * iz / l
    k = np.array([[x1, 0, 0, 0, 0, 0, -x1, 0, 0, 0, 0, 0],[0, x2, 0, 0, 0, x3, 0, -x2, 0, 0, 0, x3], \
    [0, 0, x4, 0, -x5, 0, 0, 0, -x4, 0, -x5, 0],[0, 0, 0, x6, 0, 0, 0, 0, 0, -x6, 0, 0], \
    [0, 0, -x5, 0, x7, 0, 0, 0, x5, 0, x8, 0],[0, x3, 0, 0, 0, x9, 0, -x3, 0, 0, 0, x10], \
    [-x1, 0, 0, 0, 0, 0, x1, 0, 0, 0, 0, 0],[0, -x2, 0, 0, 0, -x3, 0, x2, 0, 0, 0, -x3], \
    [0, 0, -x4, 0, x5, 0, 0, 0, x4, 0, x5, 0],[0, 0, 0, -x6, 0, 0, 0, 0, 0, x6, 0, 0], \
    [0, 0, -x5, 0, x8, 0, 0, 0, x5, 0, x7, 0],[0, x3, 0, 0, 0, x10, 0, -x3, 0, 0, 0, x9]])
    y = k

    return y