import numpy as np
from math import *

def frame3d_assembly_mat(k,kk,node_number,size_ele_mat):
    #生成整体刚度矩阵，kk为对应单元的刚度矩阵，node_number = [i;j]。
    i = node_number[0]
    j = node_number[1]
    dof = np.array([6 * i - 6,6 * i - 5,6 * i - 4,6 * i - 3,6 * i - 2,6 * i-1, \
                    6 * j - 6, 6 * j - 5, 6 * j - 4, 6 * j - 3, 6 * j - 2, 6 * j - 1])
    for m in range(size_ele_mat):
        for n in range(size_ele_mat):
            k[dof[m], dof[n]] = k[dof[m], dof[n]] + kk[m, n]
    y = k

    return y