import numpy as np
from math import *

def frame2d_ele_coordinate(node_coordinate, node_number):


    # 生成每个单元对应的节点坐标的列阵。
    eles = node_number.shape[1]
    k = np.zeros((4, eles))
    for i in range(eles):
        for j in range(2):
            k[2 * j: 2 * j + 2, i] = node_coordinate[:, node_number[j, i]-1]
    y = k

    return y
