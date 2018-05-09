import numpy as np
from math import *

def frame2d_solve_unknown(k,unknown_u_index,f,u,nodes):

    #利用分块矩阵计算未知的变量，力和位移。
    answer = []
    kkk = k.copy()   #此处变换需要用深拷贝。
    a = unknown_u_index
    n =a.shape[0]
    b = np.zeros((3 * nodes - n, 1))
    j = 0
    for i in range(3* nodes):
        if not(i+1 in a):
            b[j,0] = i+1
            j = j + 1

    #需要先后进行矩阵的行变换和列变换。
    for i in range(n):
        kkk[i,:]=k[int(a[i,0]-1),:]
    for i in range(b.shape[0]):
        kkk[n + i,:]=k[int(b[i,0]-1),:]

    kk = kkk.copy()
    for i in range(n):
        kk[:, i]=kkk[:, int(a[i,0]-1)]
    for i in range(b.shape[0]):
        kk[:, n + i]=kkk[:, int(b[i,0]-1)]

    #分块后的矩阵。
    k1 = kk[0:n, 0: n]
    k2 = kk[0:n, n : 3 * nodes]
    k3 = kk[n :3 * nodes, 0: n]
    k4 = kk[n :3 * nodes, n : 3 * nodes]

    ff=np.array((f - np.dot(k2, u)).transpose()[0],dtype=float)
    u_solved = np.linalg.solve(k1,ff)
    u_solved.shape = (u_solved.shape[0], 1)
    f_solved = np.array(np.dot(k4,u)+np.dot(k3,u_solved),dtype=float)

    #通过运算解出的位移、节点力分量。
    answer.append(u_solved)
    answer.append(f_solved)

    #F.U为总节点位移、力的列阵。
    F=np.zeros((3*nodes,1),dtype=float)
    U=np.zeros((3*nodes,1),dtype=float)
    index_f_unknown = 0
    index_u_known = 0
    index_f_known = 0
    index_u_unknown = 0
    for i in range (3 * nodes):
        #位移为已知量，力为未知量的编号对应写入分量值。
        if not(i+1 in a) :
            F[i] = f_solved[index_f_unknown,0]
            index_f_unknown = index_f_unknown + 1
            U[i] = u[index_u_known,0]
            index_u_known = index_u_known + 1
        #位移为未知量，力为已知量的编号。
        else :
            F[i] = f[index_f_known,0]
            index_f_known = index_f_known + 1
            U[i] = u_solved[index_u_unknown,0]
            index_u_unknown = index_u_unknown + 1
    answer.append(U)
    answer.append(F)

    answer.append(k1)
    answer.append(k2)
    answer.append(k3)
    answer.append(k4)

    return answer
