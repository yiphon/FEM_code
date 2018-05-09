import numpy as np
from math import *

def frame2d_assembly_mat(k, kk, node_number):

    #生成整体刚度矩阵，kk为对应单元的刚度矩阵，node_number=[i;j]。
    i=node_number[0]
    j=node_number[1]
    dof=np.array([3*i-3,3*i-2,3*i-1,3*j-3,3*j-2,3*j-1])
    for m in range(6):
        for n in range(6):
            k[dof[m],dof[n]]=k[dof[m],dof[n]]+kk[m,n]
    y=k

    return y