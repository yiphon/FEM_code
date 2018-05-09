function y = frame3d_ele_stiff_mat( r,k_ele )
%FRAME3D_ELE_STIFF_MAT Summary of this function goes here
%   Detailed explanation goes here
%   生成总体坐标系下的单元刚度矩阵，r与k_ele为转移矩阵和局部坐标系下的单元刚度矩阵。
y=r'*k_ele*r;

end

