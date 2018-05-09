import numpy as np
from math import *

def frame2d_ele_info(node_coordinate):

    #利用各单元对应节点坐标，完成ele_info矩阵。
    eles = node_coordinate.shape[1]
    a = np.zeros((3, eles))
    for i in range(eles):
        x1 = node_coordinate[0, i]
        y1 = node_coordinate[1, i]
        x2 = node_coordinate[2, i]
        y2 = node_coordinate[3, i]
        l = sqrt(pow((x2 - x1),2) + pow((y2 - y1) ,2))
        a[0, i] = l
        c = (x2 - x1) / l
        s = (y2 - y1) / l
        a[1, i] = c
        a[2, i] = s

    return a