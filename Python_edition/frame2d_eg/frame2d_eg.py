import numpy as np
from math import *
import pandas as pd
import matplotlib.pyplot as plt
from xlwt import *

from fem.frame2d_ele_coordinate import *
from fem.frame2d_ele_info import *
from fem.frame2d_r_mat import *
from fem.frame2d_k_ele_mat import *
from fem.frame2d_ele_stiff_mat import *
from fem.frame2d_assembly_mat import *
from fem.frame2d_work_equi_distributed import *
from fem.frame2d_solve_unknown import *
from fem.frame2d_beam_deflection_plot import *

def main():
    # -----------------数据的读入与有效提取。--------------------------------------------------------------------------#
    a = pd.read_excel('frame2d_eg_input.xlsx', sheet_name='Sheet1')
    b = a
    for i in a.index:
        if type(a.iloc[i, 0]) == str:
            b = b.drop([i])
    s = b.values  # s为处理后得到的源数据。之后再进行提取。

    eles = s[0, 0]  #单元个数。
    nodes = s[0, 1] #节点个数。
    dof=s[0,2]      #单个节点的自由度。
    n_for_e=int(s[0,3])  #每个单元的节点个数。
    size_ele_mat=dof*n_for_e        #单元刚度矩阵的规模。
    e = s[1, 0:eles]
    a = s[2, 0:eles]
    iz = s[3, 0:eles]
    node_number = s[4:4 + eles, 0:2].transpose()  # 每个单元对应节点的编号。
    node_coordinate = s[4 + eles:4 + eles + nodes, 0: 2].transpose()
    # node_coordinate为各单元节点对应的坐标。
    node_coordinate = frame2d_ele_coordinate(node_coordinate, node_number)

    # unknown_u_index为未知位移节点的编号。
    unknown_u_index = s[4 + eles + nodes, :]
    unknown_u_index = unknown_u_index[~pd.isnull(unknown_u_index)]
    unknown_u_index.shape = (unknown_u_index.shape[0], 1)

    # f_nodal为已知的节点力。
    f_nodal = s[5 + eles + nodes, 0:unknown_u_index.shape[0]]
    f_nodal = f_nodal[~pd.isnull(f_nodal)]
    f_nodal.shape = (f_nodal.shape[0], 1)

    # u_known为已知的节点位移。
    u_known = s[6 + eles + nodes, 0:dof * nodes - unknown_u_index.shape[0]]
    u_known = u_known[~pd.isnull(u_known)]
    u_known.shape = (u_known.shape[0], 1)

    # 以下为分布力对应的参数。
    dis_f_val = s[7 + eles + nodes, 0]  # f分布力的大小。
    dis_f_ele = s[7 + eles + nodes, 1]  # 分布力作用的单元编号。
    dis_f_range = s[7 + eles + nodes, 2]  # 分布力作用的范围。

    # 以下完成单元信息矩阵，每列依次是e\a\iz\l\cos\sin。
    ele_info = np.zeros((size_ele_mat, eles),dtype=float)
    ele_info[0, :] = e
    ele_info[1, :] = a
    ele_info[2, :] = iz
    ele_info[3: 6,:]=frame2d_ele_info(node_coordinate)



    #-----------------------单元刚度矩阵的建立。-----------------------------------------------------------------------#

    # r为单元转移矩阵的集合，每单元六行六列。
    r = np.zeros((eles,size_ele_mat, size_ele_mat))
    for i in range(eles):
        r[i,:,:]=frame2d_r_mat(ele_info[4: 6, i])


    #k_ele是单元在局部坐标系下的刚度矩阵，每单元六行六列。
    k_ele = np.zeros((eles,size_ele_mat, size_ele_mat))
    for i in range(eles):
        k_ele[i,:,:]=frame2d_k_ele_mat(ele_info[0: 4, i])

    # kk是整体坐标下单元刚度矩阵的集合，每单元六行六列。
    kk = np.zeros((eles,size_ele_mat,size_ele_mat))
    for i in range(eles):
        kk[i,:,:]=frame2d_ele_stiff_mat(r[i,:,:], k_ele[i,:,:])


    #----------------------总体方程的求解。----------------------------------------------------------------------------#
    # k是总体刚度矩阵。
    k = np.zeros((nodes * dof, nodes * dof))
    for i in range(eles):
        k = frame2d_assembly_mat(k, kk[i, :, :], node_number[:, i])

    # equi为等效节点载荷列阵。
    equi = frame2d_work_equi_distributed(dis_f_range, dis_f_val)
    f_known = f_nodal
    f_known[np.where(unknown_u_index == dof * node_number[0, dis_f_ele - 1]-1)]= f_known[np.where(unknown_u_index == dof * node_number[0, dis_f_ele - 1]-1)][0]+equi[0,0]
    f_known[np.where(unknown_u_index == dof * node_number[0, dis_f_ele - 1] )]=f_known[np.where(unknown_u_index == dof * node_number[0, dis_f_ele - 1])][0]+equi[1,0]
    f_known[np.where(unknown_u_index == dof * node_number[1, dis_f_ele - 1] - 1)]=f_known[np.where(unknown_u_index == dof * node_number[1, dis_f_ele - 1]-1)][0]+equi[2,0]
    f_known[np.where(unknown_u_index == dof * node_number[1, dis_f_ele - 1] )]=f_known[np.where(unknown_u_index == dof * node_number[1, dis_f_ele - 1])][0]+equi[3,0]

    #利用分块矩阵计算未知的变量，力和位移。
    unknown_var = frame2d_solve_unknown(k, unknown_u_index, f_known, u_known, nodes)
    ua = np.array(unknown_var[0],dtype=float)
    fc = np.array(unknown_var[1],dtype=float)
    U = np.array(unknown_var[2],dtype=float)
    F = np.array(unknown_var[3],dtype=float)
    kcc = np.array(unknown_var[4],dtype=float)
    kca = np.array(unknown_var[5],dtype=float)
    kac = np.array(unknown_var[6],dtype=float)
    kaa = np.array(unknown_var[7],dtype=float)


    #---------------------------后处理。-------------------------------------------------------------------------------#

    #u为每单元两节点的位移阵列。
    u = np.zeros((size_ele_mat, eles),dtype=float)
    for i in range(eles):
        u[0:3, i]=np.array(U[node_number[0, i] * 3 -3: node_number[0, i] * 3,0],dtype=float)
        u[3:6,i]=np.array(U[node_number[1, i] * 3 - 3: node_number[1, i] * 3,0],dtype=float)

    #f为每单元两节点的受力阵列，按以下对每个单元操作的步骤后，还应减去等效节点载荷。
    f_nodal = np.zeros((size_ele_mat, eles),dtype=float)
    for i in range(eles):
        f_nodal[:, i]=np.dot(np.dot(k_ele[i,:,:],r[i,:,:]),u[:, i])
        if i == dis_f_ele-1:
            equi = np.array([0,equi[0,0],equi[1,0],0,equi[2,0],equi[3,0]],dtype=float)
            f_nodal[:, i]=f_nodal[:, i]-equi

    #绘出各单元挠度转角图像。
    fig1,axes=plt.subplots(eles,2)
    for i in range(eles):
        frame2d_beam_deflection_plot(ele_info[3, i], u[[1, 2, 4, 5], i], eles, i,axes[i,0],axes[i,1])
    # plt.show()


    #=============文件写入。============================================================================================

    size_ele_mat
    wb=Workbook(encoding='utf-8')
    sh=wb.add_sheet('frame_2d')

    sh.write(0,0,'element stiffness matrices:')
    for i in range(eles):
        for j in range(size_ele_mat):
            for m in range(size_ele_mat):
                sh.write(i*7+j+1,m,kk[i,j,m])

    sh.write((size_ele_mat+1)*eles+1,0,'global stiffness matrix:')
    for j in range(dof * nodes):
        for m in range(dof * nodes):
            sh.write((size_ele_mat+1)*eles+2+j,m,k[j, m])

    sh.write((size_ele_mat+1)*eles+dof*nodes+3,0, 'K_cc=')
    for i in range(kcc.shape[0]):
        for j in range(kcc.shape[1]):
            sh.write((size_ele_mat+1)*eles+dof*nodes+4+i,j,kcc[i, j])

    sh.write((size_ele_mat+1)*eles+dof*nodes+5+kcc.shape[0],0, 'K_ac=')
    for i in range(kac.shape[0]):
        for j in range(kac.shape[1]):
            sh.write(7*eles+3*nodes+6+kcc.shape[0]+i,j,kac[i, j])

    sh.write((size_ele_mat+1)*eles+2*dof*nodes+7,0,'K_ca=')
    for i in range(kca.shape[0]):
        for j in range(kca.shape[1]):
            sh.write((size_ele_mat+1)*eles+2*dof*nodes+8+i,j,kca[i, j])

    sh.write((size_ele_mat+1)*eles+2*dof*nodes+9+kca.shape[0],0, 'K_aa=')
    for i in range(kaa.shape[0]):
        for j in range(kaa.shape[1]):
            sh.write((size_ele_mat+1)*eles+2*dof*nodes+10+kca.shape[0]+i,j,kaa[i, j])

    sh.write((size_ele_mat+1)*eles+3*dof*nodes+11,0, 'U_a=')
    for i in range(ua.shape[0]):
        sh.write((size_ele_mat+1)*eles+3*dof*nodes+12,i,ua[i,0])

    sh.write((size_ele_mat+1)*eles+3*dof*nodes+14,0, 'F_c=')
    for i in range(fc.shape[0]):
        sh.write((size_ele_mat+1)*eles+3*dof*nodes+15,i,fc[i,0])

    sh.write((size_ele_mat+1)*eles+3*dof*nodes+17,0, 'U=')
    for i in range(U.shape[0]):
        sh.write((size_ele_mat+1)*eles+3*dof*nodes+18,i,U[i,0])

    sh.write((size_ele_mat+1)*eles+3*dof*nodes+20,0, 'F=')
    for i in range(U.shape[0]):
        sh.write((size_ele_mat+1)*eles+3*dof*nodes+21,i,F[i,0])

    sh.write((size_ele_mat+1)*eles+3*dof*nodes+23,0, 'f_nodal=')
    for i in range(f_nodal.shape[0]):
        for j in range(f_nodal.shape[1]):
            sh.write((size_ele_mat+1)*eles+3*dof*nodes+24+i,2*j,f_nodal[i, j])

    wb.save('frame2d_eg_output_py.xls')

#=======================================================================================================================

if __name__ =='__main__':
    main()
