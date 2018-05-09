import numpy as np
from math import *

def frame2d_r_mat(ele_info):

    #生成单元转换矩阵。
    c = ele_info[0]
    s = ele_info[1]
    rr = np.array([[c, s, 0],[-s, c, 0],[0, 0, 1]])
    r1=np.column_stack((rr,np.zeros((3, 3))))
    r2=np.column_stack((np.zeros((3, 3)),rr))
    r=np.row_stack((r1,r2))
    y = r

    return y