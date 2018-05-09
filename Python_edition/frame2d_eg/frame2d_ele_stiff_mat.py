import numpy as np
from math import *

def frame2d_ele_stiff_mat(r,k_ele):

    #完成整体坐标系下单元刚度矩阵，r与k是对应单元的转移矩阵和局部坐标单元刚度矩阵。
    a=np.dot(r.transpose(),k_ele)
    a=np.dot(a,r)
    y = a

    return y