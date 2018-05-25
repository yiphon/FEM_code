import numpy as np
from math import *


def frame3d_ele_info(node_coordinate):
    # 利用每个单元的节点坐标列阵，把单元长度、正余弦补充完整。
    eles = node_coordinate.shape[1]
    a = np.zeros((4, eles))
    for i in range(eles):
        x1 = node_coordinate[0, i]
        y1 = node_coordinate[1, i]
        z1 = node_coordinate[2, i]
        x2 = node_coordinate[3, i]
        y2 = node_coordinate[4, i]
        z2 = node_coordinate[5, i]
        l = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
        a[0, i] = l
        cxx = (x2 - x1) / l
        cyx = (y2 - y1) / l
        czx = (z2 - z1) / l
        a[1, i] = cxx
        a[2, i] = cyx
        a[3, i] = czx
    y = a

    return y
